#![allow(dead_code)]
//! Canonical chemical equilibrium solver using log-mole formulation and reaction extent approach.
//!
//! # Purpose
//!
//! This module is the **central orchestrator** of the chemical equilibrium calculation
//! pipeline. It owns all equilibrium data (elemental composition, stoichiometry, Gibbs
//! functions, solver settings), manages the solve lifecycle, and publishes results.
//! The module solves chemical equilibrium problems by **minimizing Gibbs free energy**
//! using reaction extents in log-mole space.
//!
//! # Physical and Mathematical Background
//!
//! ## Log-Mole Formulation
//!
//! The solver uses **logarithmic mole variables** `y_i = ln(n_i)` instead of physical
//! moles `n_i`. This transformation:
//!
//! - Automatically enforces `n_i > 0` (physical constraint).
//! - Handles trace species with very small concentrations without underflow.
//! - Linearizes the mass-action equations near equilibrium.
//!
//! The residual vector `F(y)` has two blocks:
//!
//! 1. **Element conservation** (E equations, one per element):
//!    ```text
//!    F_elem = A^T · n(y) - b_0
//!    ```
//!    where `A` is the `(species × elements)` composition matrix, `n(y) = exp(y)`
//!    are physical moles, and `b_0` are the total element moles from the initial state.
//!
//! 2. **Reaction affinity** (R equations, one per independent reaction):
//!    ```text
//!    F_rxn = ν^T · (μ(y) / (R·T))
//!    ```
//!    where `ν` is the `(species × reactions)` stoichiometric matrix from SVD,
//!    and `μ(y) = μ°(T) + R·T·ln(a_i(y))` is the chemical potential.
//!
//! ## SVD Reaction Basis
//!
//! Independent reactions are computed via **SVD of the element composition matrix**:
//!
//! ```text
//! A = U · Σ · V^T
//! ```
//!
//! The nullspace of `A^T` (columns of `V` corresponding to zero singular values)
//! forms the stoichiometric matrix `ν`. This guarantees that reaction extents
//! automatically conserve elements.
//!
//! ## Scaling
//!
//! Row scaling normalizes residual rows by their characteristic magnitude,
//! improving numerical conditioning for systems with species at very different
//! concentration levels.
//!
//! # Architecture and Dataflow
//!
//! ```text
//!   ┌─────────────────────────────────────────────────────────────────────┐
//!   │                        EquilibriumLogMoles                          │
//!   │                                                                     │
//!   │  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────────┐  │
//!   │  │  Problem Data │  │  Solver Data │  │  Publication State       │  │
//!   │  │  ─ elem_comp  │  │  ─ settings  │  │  ─ solution (log_moles)  │  │
//!   │  │  ─ stoich_mat │  │  ─ params    │  │  ─ moles (physical)     │  │
//!   │  │  ─ gibbs      │  │  ─ policy    │  │  ─ solve_report         │  │
//!   │  │  ─ phases     │  │  ─ budget    │  │  ─ keq_validation       │  │
//!   │  │  ─ n0         │  │  ─ seed      │  │  ─ mole_table           │  │
//!   │  │  ─ T, P, p0   │  │              │  │                          │  │
//!   │  └──────┬───────┘  └──────┬───────┘  └──────────┬───────────────┘  │
//!   └─────────┼─────────────────┼──────────────────────┼──────────────────┘
//!             │                 │                      │
//!             v                 v                      v
//!   ┌─────────────────────────────────────────────────────────────────────┐
//!   │  Solve Pipeline                                                     │
//!   │                                                                     │
//!   │  1. check_task() ──> validate dimensions and settings               │
//!   │  2. stage_equilibrium_system() ──> build SVD reaction basis         │
//!   │  3. build_solve_contract() ──> create residual/Jacobian closures    │
//!   │  4. solver_impl() ──> run backend cascade                           │
//!   │     ├── solve_backend_cascade() ──> try backends in order           │
//!   │     │   ├── Legacy(LM/NR/TR) ──> equilibrium_legacy_backend         │
//!   │     │   └── RustedSciThe ──> equilibrium_rst_backend                │
//!   │     └── validate_equilibrium_candidate() ──> acceptance gate        │
//!   │  5. publish_solve_candidate() ──> store accepted solution           │
//!   │  6. run_keq_cross_validation() ──> independent K_eq check           │
//!   └─────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! ## Temperature Range Pipeline
//!
//! ```text
//!   solve_for_T_range(T_start, T_end, T_step)
//!     │
//!     ├── build_temperature_ranges() ──> split range into segments
//!     ├── build_temperature_gibbs_cache() ──> precompute Gibbs functions
//!     ├── build_temperature_point_index() ──> map T → range segment
//!     │
//!     └── for each temperature point:
//!           ├── continuation_seed_for_point() ──> seed from previous solution
//!           ├── solve_temperature_point_from_seed() ──> solve at this T
//!           └── publish_solve_candidate() ──> store result
//! ```
//!
//! # Key Structures
//!
//! | Structure | Role |
//! |-----------|------|
//! | [`EquilibriumLogMoles`] | Central orchestrator — owns all data and methods |
//! | [`SolverParams`] | Numerical solver tuning parameters (λ, tol, max_iter, etc.) |
//! | [`EquilibriumSolverSettings`] | High-level settings (backend policy, scaling, K_eq validation) |
//! | [`EquilibriumSolveContract`] | Bundles residual, Jacobian, and feasibility closures |
//! | [`EquilibriumSolveCandidate`] | One accepted candidate (solution + reports) |
//! | [`TemperatureWorkerSeed`] | Snapshot of solver state for temperature-range parallelism |
//! | [`TemperatureSolveSnapshot`] | Result of one temperature point solve |
//! | [`TemperatureSolveFailure`] | Error information for a failed temperature point |
//! | [`Phase`] | Phase descriptor (activity model + species indices) |
//! | [`Solvers`] | Enum: LM, NR, TR |
//!
//! # Key Free Functions
//!
//! | Function | Purpose |
//! |----------|---------|
//! | [`compute_species_moles`] | Converts log-moles `y` to physical moles `n = exp(y)` |
//! | [`compute_element_totals`] | Computes `A^T · n` — element totals from composition matrix |
//! | [`reaction_standard_gibbs`] | Computes `ν^T · g(T)` — standard Gibbs for each reaction |
//! | [`equilibrium_scaling`] | Computes row scaling factors for residual/Jacobian |
//! | [`species_to_phase_map`] | Maps each species index to its phase index |
//! | [`reaction_phase_stoichiometry`] | Aggregates stoichiometry by phase |
//! | [`evaluate_equilibrium_logmole_residual`] | Evaluates the full residual vector |
//! | [`evaluate_equilibrium_logmole_jacobian`] | Evaluates the full Jacobian matrix |
//! | [`scaled_residual`] / [`scaled_jacobian`] | Wraps residual/Jacobian with row scaling |
//! | [`validate_logmole_system_dimensions`] | Validates all matrix/vector dimensions |
//! | [`validate_residual_conditions`] | Validates T, P, p0 for residual evaluation |
//!
//! # Examples
//!
//! ```rust, ignore
//! use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::*;
//!
//! // Create equilibrium solver for gas-phase system
//! let mut solver = gas_solver(
//!     vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
//!     1000.0, 101325.0, Solvers::LM, Some("info"), true
//! );
//! solver.solve().unwrap();
//! let solution = solver.accepted_solution().unwrap();
//! ```
//!
//! # Non-obvious Details
//!
//! - **Log-mole variables**: `y_i = ln(n_i)`. The solver works in log-mole space;
//!   `solver.solution` contains log-moles, not physical moles. Use
//!   [`compute_species_moles`] to convert.
//! - **SVD reaction basis**: The stoichiometric matrix is the nullspace of `A^T`,
//!   computed via SVD. This guarantees element conservation by construction.
//! - **Row scaling**: [`equilibrium_scaling`] computes per-row scale factors based on
//!   element totals and reaction standard Gibbs. Scaling is optional but recommended
//!   for systems with large concentration disparities.
//! - **Temperature continuation**: When solving over a temperature range, the solution
//!   at each point seeds the next point (unless `ContinuationSeedPolicy::IndependentPerPoint`
//!   is set). This dramatically improves convergence.
//! - **Backend cascade**: If one solver backend fails, the cascade tries the next
//!   backend in order, subject to the iteration budget.
//!
//! # Related Modules
//!
//! - [`equilibrium_nonlinear`](super::equilibrium_nonlinear) — numerical solver implementations
//! - [`equilibrium_workflows`](super::equilibrium_workflows) — phase control and convenience workflows
//! - [`equilibrium_problem`](super::equilibrium_problem) — typed input boundary
//! - [`equilibrium_validation`](super::equilibrium_validation) — backend-independent acceptance gate
//! - [`equilibrium_solver_policy`](super::equilibrium_solver_policy) — backend selection and cascade budget
//! - [`equilibrium_rst_backend`](super::equilibrium_rst_backend) — RustedSciThe symbolic backend
//! - [`equilibrium_legacy_backend`](super::equilibrium_legacy_backend) — legacy LM/NR/TR adapter
//! - [`equilibrium_activity`](super::equilibrium_activity) — activity models
//! - [`equilibrium_ids`](super::equilibrium_ids) — typed index wrappers
//! - [`equilibrium_constant_cross_validation`](super::equilibrium_constant_cross_validation) — K_eq cross-validation
//! - [`equilibrium_temperature_postprocessing`](super::equilibrium_temperature_postprocessing) — postprocessing
//!
use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::{
    PhaseActivityModel, phase_activity_models,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::{
    BackendSolveRequest, EquilibriumNonlinearBackend,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_cross_validation::{
    EquilibriumConstantCrossValidationStatus, EquilibriumConstantCrossValidationTolerances,
    classify_equilibrium_constant_cross_validation,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::EquilibriumConstantProblem;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::{
    EquilibriumConstantSolver, EquilibriumConstantSolverMode,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationMode, EquilibriumConstantValidationTolerances,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionBasis, ReactionExtentError, compute_reaction_basis,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    DEFAULT_TRACE_MOLE_FLOOR, EquilibriumConditions, EquilibriumProblem, EquilibriumSolution,
    LogMolesInitialGuess, PreparedEquilibriumProblem, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, prepare_rst_symbolic_problem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    EquilibriumSolveReport, SolverAttemptFailureKind, SolverAttemptOutcome, SolverAttemptReport,
    SolverBackend, SolverCascadeBudget, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    TemperaturePostprocessingPolicy, TemperaturePostprocessingResult,
    postprocess_temperature_series,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::{
    EquilibriumAcceptanceCriteria, EquilibriumCandidateReport, EquilibriumCandidateResiduals,
    validate_equilibrium_candidate,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
    PhaseControlledSolveReport, PhaseManager, multiphase_equilibrium_residual_generator_sym,
};
use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::SubsDataError;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_functions::Jacobian;
use log::{error, info, warn};
use nalgebra::{DMatrix, DVector};
use prettytable::{Cell, Row, Table};
use std::collections::{HashMap, HashSet};
use std::default::Default;
use std::f64;
use std::rc::Rc;
use std::time::Instant;
/// Universal gas constant in J/(mol·K)
pub const R: f64 = 8.314;

/// Configuration parameters for numerical solvers.
///
/// These parameters control the convergence behavior of the nonlinear solvers
/// (LM, NR, Trust Region). They are validated by
/// [`EquilibriumSolverSettings::validate()`] before any solver is constructed.
#[derive(Clone)]
pub struct SolverParams {
    /// Maximum number of nonlinear iterations before giving up.
    pub max_iter: usize,
    /// Convergence tolerance for the residual L2 norm.
    /// The solver stops when `||F(x)|| < tol`.
    pub tol: f64,
    /// Initial damping parameter for the Levenberg-Marquardt solver.
    /// Controls the blend between Gauss-Newton (λ → 0) and gradient descent (λ → ∞).
    pub lambda: f64,
    /// Minimum step size in line search before declaring failure.
    /// If the Armijo condition cannot be satisfied with α ≥ alpha_min, the solver stops.
    pub alpha_min: f64,
    /// Initial trust-region radius for the Trust Region solver.
    pub delta_init: f64,
    /// Maximum trust-region radius. Must be ≥ delta_init.
    pub delta_max: f64,
    /// Acceptance threshold for the trust-region ratio ρ.
    /// If ρ > eta, the step is accepted; if ρ < eta/4, the radius shrinks.
    /// Must lie in the open interval (0, 1).
    pub eta: f64,
}

/// Backend-specific solver controls separated from the physical equilibrium data.
///
/// The `EquilibriumLogMoles` facade still carries domain data, sweep state, and
/// published results. This structure groups the numerical policy knobs that are
/// specific to backend selection and continuation behavior.
#[derive(Clone)]
pub struct EquilibriumSolverSettings {
    /// Numerical solver parameters (max iterations, tolerances, damping, trust-region radii).
    pub solver_params: SolverParams,
    /// Preferred legacy solver backend (LM, NR, or TR) used when `solver_policy` is `None`.
    pub solver: Solvers,
    /// Whether to enable row scaling of the residual and Jacobian.
    /// When `true`, each residual/Jacobian row is divided by its scale factor.
    pub scaling_flag: bool,
    /// Typed policy for generating implicit log-mole seeds from physical
    /// initial moles when the caller did not provide an explicit seed.
    pub trace_seed_policy: TraceSpeciesSeedPolicy,
    /// Optional explicit backend policy. `None` retains the historical order
    /// beginning with [`Self::solver`] while callers migrate to typed policies.
    pub solver_policy: Option<SolverPolicy>,
    /// Optional resource limits for a solver cascade. `None` derives a budget
    /// from the current policy and legacy per-backend iteration setting.
    pub solver_budget: Option<SolverCascadeBudget>,
    /// Controls whether the accepted solution is additionally checked by the
    /// independent equilibrium-constant validator.
    pub keq_validation_mode: EquilibriumConstantValidationMode,
    /// Numerical tolerances for the independent cross-validation report.
    pub keq_validation_tolerances: EquilibriumConstantCrossValidationTolerances,
    /// Explicit warm-start behavior for sequential and chunked temperature
    /// sweeps. The fully parallel `par2` sweep is always independent.
    pub continuation_seed_policy: ContinuationSeedPolicy,
}

impl Default for EquilibriumSolverSettings {
    fn default() -> Self {
        Self {
            solver_params: SolverParams::default(),
            solver: Solvers::LM,
            scaling_flag: false,
            trace_seed_policy: TraceSpeciesSeedPolicy::Absolute {
                floor: DEFAULT_TRACE_MOLE_FLOOR,
            },
            solver_policy: None,
            solver_budget: None,
            keq_validation_mode: EquilibriumConstantValidationMode::Off,
            keq_validation_tolerances: EquilibriumConstantCrossValidationTolerances::default(),
            continuation_seed_policy: ContinuationSeedPolicy::default(),
        }
    }
}

impl EquilibriumSolverSettings {
    /// Validates numerical controls before any residual, Jacobian, or backend
    /// is constructed.
    ///
    /// Keeping this check with the settings bundle makes the canonical
    /// one-shot API and the temporary mutable facade reject the same malformed
    /// solver request at the same boundary.
    pub fn validate(&self) -> Result<(), ReactionExtentError> {
        let params = &self.solver_params;
        if params.max_iter == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_params.max_iter",
                message: "max_iter must be greater than zero".to_string(),
            });
        }
        for (field, value) in [
            ("solver_params.tol", params.tol),
            ("solver_params.lambda", params.lambda),
            ("solver_params.alpha_min", params.alpha_min),
            ("solver_params.delta_init", params.delta_init),
            ("solver_params.delta_max", params.delta_max),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field,
                    message: "value must be finite and strictly positive".to_string(),
                });
            }
        }
        if params.delta_max < params.delta_init {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_params.delta_max",
                message: "delta_max must be greater than or equal to delta_init".to_string(),
            });
        }
        if !params.eta.is_finite() || !(0.0..1.0).contains(&params.eta) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_params.eta",
                message: "eta must be finite and lie in the open interval (0, 1)".to_string(),
            });
        }

        if let Some(policy) = &self.solver_policy {
            if policy.ordered_backends().is_empty() {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "solver_policy",
                    message: "a solver policy must contain at least one backend".to_string(),
                });
            }
        }
        if let Some(budget) = self.solver_budget {
            if budget.max_attempts == 0
                || budget.max_iterations_per_attempt == 0
                || budget.max_total_iterations == 0
            {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "solver_budget",
                    message: "attempt, per-attempt, and total iteration limits must be positive"
                        .to_string(),
                });
            }
        }
        self.keq_validation_tolerances.validate()?;
        Ok(())
    }
}

/// Available numerical solvers for equilibrium calculations.
///
/// These are the legacy hand-written solvers. New code should prefer
/// [`SolverBackend::RustedSciThe`](crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverBackend::RustedSciThe)
/// via the [`SolverPolicy`](crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy) system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Solvers {
    /// Levenberg-Marquardt solver with adaptive damping.
    LM,
    /// Newton-Raphson solver with line search and feasibility constraints.
    NR,
    /// Trust Region solver with adaptive radius (dogleg method).
    TR,
}

/// Defines how a temperature sweep obtains the initial log-mole iterate for
/// each requested temperature point.
///
/// `PreviousAccepted` is a continuation method: an accepted point warms the
/// next point. `IndependentPerPoint` always reuses the seed configured before
/// the sweep, which is useful for reproducible point-by-point comparisons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuationSeedPolicy {
    PreviousAccepted,
    IndependentPerPoint,
}

impl Default for ContinuationSeedPolicy {
    fn default() -> Self {
        Self::PreviousAccepted
    }
}

impl Default for SolverParams {
    fn default() -> Self {
        Self {
            max_iter: 50,
            tol: 1e-6,
            lambda: 1e-3,
            alpha_min: 1e-6,
            delta_init: 1.0,
            delta_max: 100.0,
            eta: 0.1,
        }
    }
}
/// Central orchestrator for chemical equilibrium calculations.
///
/// This struct owns all equilibrium data (elemental composition, stoichiometry,
/// Gibbs functions, solver settings), manages the solve lifecycle, and publishes
/// results. It is the primary entry point for the canonical equilibrium solver.
///
/// # Dataflow
///
/// See the [module-level documentation](self) for a complete dataflow diagram.
///
/// # Warning
///
/// This struct has **22 public fields** (see SourceCraft Diagnostics item A.4
/// in [`TODO_ANALYSIS.md`](https://todo)). Many fields are interdependent and
/// should be modified through methods rather than directly. Prefer using
/// [`EquilibriumProblem`](super::equilibrium_problem::EquilibriumProblem) for
/// constructing new problems and [`from_problem()`](Self::from_problem) for
/// loading them into the solver.
pub struct EquilibriumLogMoles {
    /// Substance database with thermochemical data for all species.
    pub subs_data: SubsData,
    /// Element composition matrix `A` of shape `(species × elements)`.
    /// `A[i, j]` = number of atoms of element `j` in species `i`.
    pub elem_composition: DMatrix<f64>,
    /// SVD-derived reaction basis containing the stoichiometric nullspace.
    pub reaction_basis: ReactionBasis,
    /// Optional explicit initial guess in log-mole space (`y_i = ln(n_i)`).
    /// If `None`, the solver derives a seed from `n0` and the trace policy.
    pub initial_guess: Option<Vec<f64>>,
    /// System pressure in Pascals (Pa).
    pub P: f64,
    /// Temperature in Kelvin (K).
    pub T: f64,
    /// Stoichiometric matrix `ν` of shape `(species × reactions)`.
    /// Columns are independent reaction directions from the SVD nullspace.
    pub stoich_matrix: DMatrix<f64>,
    /// Initial physical moles `n_i^0` for each species.
    pub n0: Vec<f64>,
    /// Standard Gibbs free energy functions `g_i(T)` for each species, in J/mol.
    pub gibbs: Vec<GibbsFn>,
    /// Symbolic expressions for Gibbs free energy (used by RustedSciThe backend).
    pub gibbs_sym: Vec<Expr>,
    /// Phase descriptors: activity model and species indices for each phase.
    pub phases: Vec<Phase>,
    /// Maps each species index to its phase index: `species_phase[i]` = phase of species i.
    pub species_phase: Vec<usize>,
    /// Persistent phase-control mask for the current equilibrium problem.
    /// It survives outer-loop restarts and is reset only when a new problem
    /// or phase inventory is staged.
    pub(crate) phase_active_mask: Vec<bool>,
    /// Legacy temperature-only failure list retained for callers that only
    /// need a quick visual summary. Prefer [`Self::temperature_failures`].
    pub list_of_failed_T: Vec<f64>,
    /// Typed diagnostics for every temperature point whose solve was not
    /// accepted during the most recent sweep.
    pub temperature_failures: Vec<TemperatureSolveFailure>,
    /// Immutable evidence for every accepted temperature point from the most
    /// recent sweep, in deterministic temperature order.
    pub temperature_solutions: Vec<TemperatureSolveSnapshot>,
    /// Accepted solution in log-mole space: `y_i = ln(n_i)`.
    /// Published by [`publish_solve_candidate()`](Self::publish_solve_candidate).
    pub solution: Vec<f64>,
    /// Accepted solution in physical moles: `n_i = exp(y_i)`.
    /// Published simultaneously with [`solution`](Self::solution).
    pub moles: Vec<f64>,
    /// Validation evidence for the currently published solution.
    pub last_validation_report: Option<EquilibriumCandidateReport>,
    /// Backend-selection trace for the currently published solution.
    pub last_solve_report: Option<EquilibriumSolveReport>,
    /// Optional independent K_eq cross-validation status for the last accepted solve.
    pub last_keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
    /// Active-set transition evidence for the last accepted phase-controlled solve.
    pub last_phase_control_report: Option<PhaseControlledSolveReport>,
    /// Element totals vector `b_0 = A^T · n0` — total moles of each element.
    /// Used in the element conservation residual block.
    pub elements_vector: Vec<f64>,
    /// Optional log level for diagnostic output (e.g., "info", "debug").
    pub loglevel: Option<String>,
    /// Species-level epsilon for phase detection thresholds.
    pub species_eps: f64,
    /// Substate-level epsilon for phase detection thresholds.
    pub substate_eps: f64,
    /// Map from substance name to vector of moles across temperature points.
    /// Populated by [`map_of_moles_for_each_substance()`](Self::map_of_moles_for_each_substance).
    pub map_of_moles_for_each_substance: HashMap<String, Vec<f64>>,
    /// Raw temperature-range results: `(temperature, moles_at_this_T)` pairs.
    pub moles_for_T_range: Vec<(f64, Vec<f64>)>,
    /// Solver configuration: backend policy, scaling, K_eq validation, etc.
    pub solver_settings: EquilibriumSolverSettings,
    /// Phase manager for multiphase equilibrium control (hysteresis, thresholds).
    pub phase_manager: PhaseManager,
    /// Reference pressure in Pascals (typically 101325 Pa = 1 atm).
    pub p0: f64,
}

/// Immutable bundle of closures and symbolic state for one nonlinear solve.
///
/// The builder centralizes the formulas for single solves and every sweep
/// variant so backend selection cannot drift between entry points.
#[allow(dead_code)]
pub(crate) struct EquilibriumSolveContract {
    /// Unscaled residual closure: `F_raw(y) = [A^T·n(y) - b0; ν^T·μ/(R·T)]`.
    f_raw: Rc<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>,
    /// Scaled residual closure (if scaling is enabled), or same as `f_raw`.
    f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>,
    /// Optional analytical Jacobian closure `J(y) = dF/dy`.
    /// `None` for backends that compute their own Jacobian (e.g., RustedSciThe).
    j: Option<Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>>,
    /// Optional precomputed symbolic problem for the RustedSciThe backend.
    rst_problem: Option<RstPreparedProblem>,
}

/// Fully validated result of one nonlinear solve before it is published.
///
/// Phase-control restarts use this bundle transactionally: an intermediate
/// candidate may guide the next active-set decision, but it cannot overwrite
/// the last accepted public state until the outer loop reaches a fixed point.
pub(crate) struct EquilibriumSolveCandidate {
    /// Accepted solution in log-mole space: `y_i = ln(n_i)`.
    pub(crate) log_moles: Vec<f64>,
    /// Reconstructed physical moles: `n_i = exp(y_i)`.
    pub(crate) moles: Vec<f64>,
    /// Map from substance name to moles (for display/export).
    pub(crate) mole_table: HashMap<String, Vec<f64>>,
    /// Backend-independent acceptance evidence.
    pub(crate) validation_report: EquilibriumCandidateReport,
    /// Backend-selection trace for this candidate.
    pub(crate) solve_report: EquilibriumSolveReport,
    /// Optional K_eq cross-validation status.
    pub(crate) keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
}

#[allow(dead_code)]
impl EquilibriumSolveContract {
    fn residual(&self, values: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        (self.f.as_ref())(values)
    }

    fn raw_residual(&self, values: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
        (self.f_raw.as_ref())(values)
    }

    fn jacobian(&self) -> Option<&dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>> {
        self.j.as_deref()
    }
}

/// Snapshot of the state copied into a temperature-sweep worker.
///
/// This keeps the parallel sweep closures free of the main solver object, so
/// rayon never sees the `Rc`-backed fields that make the facade itself
/// non-`Sync`.
/// Snapshot of the state copied into a temperature-sweep worker.
///
/// This keeps the parallel sweep closures free of the main solver object, so
/// rayon never sees the `Rc`-backed fields that make the facade itself
/// non-`Sync`.
#[derive(Clone)]
pub(crate) struct TemperatureWorkerSeed {
    /// Cloned substance database for the worker.
    subs_data: SubsData,
    /// Element composition matrix `A` (species × elements).
    elem_composition: DMatrix<f64>,
    /// SVD-derived reaction basis.
    reaction_basis: ReactionBasis,
    /// Stoichiometric matrix `ν` (species × reactions).
    stoich_matrix: DMatrix<f64>,
    /// Initial physical moles `n_i^0`.
    n0: Vec<f64>,
    /// Phase descriptors.
    phases: Vec<Phase>,
    /// Species-to-phase mapping.
    species_phase: Vec<usize>,
    /// Element totals vector `b_0`.
    elements_vector: Vec<f64>,
    /// Optional log level for diagnostic output.
    loglevel: Option<String>,
    /// Species-level epsilon for phase detection.
    species_eps: f64,
    /// Substate-level epsilon for phase detection.
    substate_eps: f64,
    /// Solver configuration (backend policy, scaling, etc.).
    solver_settings: EquilibriumSolverSettings,
    /// System pressure in Pa.
    P: f64,
    /// Reference pressure in Pa.
    p0: f64,
}

impl TemperatureWorkerSeed {
    /// Copies the frozen seed state into a local solver instance for one temperature point.
    ///
    /// This is the core cloning mechanism for parallel temperature sweeps. It copies
    /// all problem data (element composition, reaction basis, stoichiometry, phases,
    /// solver settings) into the local solver, sets the temperature, and clears any
    /// previously published state. The worker can then call `local_solver.solve()`
    /// independently without affecting other workers or the main solver.
    pub(crate) fn apply(
        &self,
        local_solver: &mut EquilibriumLogMoles,
        temperature: f64,
        gibbs: Vec<GibbsFn>,
        initial_guess: Option<Vec<f64>>,
    ) {
        local_solver.subs_data = self.subs_data.clone();
        local_solver.elem_composition = self.elem_composition.clone();
        local_solver.reaction_basis = self.reaction_basis.clone();
        local_solver.stoich_matrix = self.stoich_matrix.clone();
        local_solver.n0 = self.n0.clone();
        local_solver.gibbs = gibbs;
        local_solver.gibbs_sym = Vec::new();
        local_solver.phases = self.phases.clone();
        local_solver.species_phase = self.species_phase.clone();
        local_solver.list_of_failed_T.clear();
        local_solver.temperature_failures.clear();
        local_solver.temperature_solutions.clear();
        local_solver.solution.clear();
        local_solver.moles.clear();
        local_solver.last_validation_report = None;
        local_solver.last_solve_report = None;
        local_solver.last_keq_validation_status = None;
        local_solver.last_phase_control_report = None;
        local_solver.elements_vector = self.elements_vector.clone();
        local_solver.loglevel = self.loglevel.clone();
        local_solver.species_eps = self.species_eps;
        local_solver.substate_eps = self.substate_eps;
        local_solver.moles_for_T_range.clear();
        local_solver.solver_settings = self.solver_settings.clone();
        local_solver.P = self.P;
        local_solver.T = temperature;
        local_solver.p0 = self.p0;
        local_solver.initial_guess = initial_guess;
    }

    /// Returns the number of species (length of `n0`).
    pub(crate) fn species_count(&self) -> usize {
        self.n0.len()
    }

    /// Returns a reference to the initial moles vector.
    pub(crate) fn n0(&self) -> &[f64] {
        &self.n0
    }

    /// Returns the pressure in Pa.
    pub(crate) fn pressure(&self) -> f64 {
        self.P
    }

    /// Returns the reference pressure in Pa.
    pub(crate) fn reference_pressure(&self) -> f64 {
        self.p0
    }

    /// Returns a reference to the solver settings.
    pub(crate) fn solver_settings(&self) -> &EquilibriumSolverSettings {
        &self.solver_settings
    }

    /// Returns a reference to the element composition matrix.
    pub(crate) fn elem_composition(&self) -> &DMatrix<f64> {
        &self.elem_composition
    }

    /// Returns a reference to the stoichiometry matrix.
    pub(crate) fn stoich_matrix(&self) -> &DMatrix<f64> {
        &self.stoich_matrix
    }

    /// Returns a reference to the phases vector.
    pub(crate) fn phases(&self) -> &[Phase] {
        &self.phases
    }

    /// Returns a reference to the species-to-phase map.
    pub(crate) fn species_phase(&self) -> &[usize] {
        &self.species_phase
    }

    /// Returns a reference to the element totals vector.
    pub(crate) fn elements_vector(&self) -> &[f64] {
        &self.elements_vector
    }

    /// Returns the species epsilon (phase detection threshold).
    pub(crate) fn species_eps(&self) -> f64 {
        self.species_eps
    }

    /// Returns the substate epsilon.
    pub(crate) fn substate_eps(&self) -> f64 {
        self.substate_eps
    }
}

/// One non-accepted temperature point from a continuation or parallel sweep.
///
/// Failure diagnostics are kept separate from successful rows so consumers
/// cannot mistake a stale solution for a result at the failed temperature.
#[derive(Debug, Clone)]
pub struct TemperatureSolveFailure {
    /// Temperature requested for this point, in K.
    pub temperature: f64,
    /// Stable diagnostic text for non-cascade and setup failures.
    pub message: String,
    /// Ordered backend trace when the failure exhausted a solver cascade.
    pub attempts: Vec<SolverAttemptReport>,
}

/// One accepted temperature point together with the evidence that admitted it.
#[derive(Debug, Clone)]
pub struct TemperatureSolveSnapshot {
    /// Temperature requested for this point, in K.
    pub temperature: f64,
    /// Accepted nonlinear coordinates `ln(n_i)` in canonical species order.
    pub log_moles: Vec<f64>,
    /// Reconstructed physical moles in the same species order.
    pub moles: Vec<f64>,
    /// Backend-independent acceptance evidence.
    pub validation: EquilibriumCandidateReport,
    /// Ordered numerical backend trace for this accepted point.
    pub solve_report: EquilibriumSolveReport,
}

/// Constructs a typed failure record from a temperature point that could not be solved.
///
/// Extracts the ordered backend cascade trace from `AllBackendsFailed` errors,
/// or creates a simple diagnostic message for non-cascade failures (e.g., invalid
/// temperature range, dimension mismatch).
pub(crate) fn temperature_failure(temperature: f64, error: &ReactionExtentError) -> TemperatureSolveFailure {
    let attempts = match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts.clone(),
        _ => Vec::new(),
    };
    TemperatureSolveFailure {
        temperature,
        message: error.to_string(),
        attempts,
    }
}

/// Selects the seed for one temperature point without mutating solver state.
///
/// Keeping this selection pure makes the continuation contract testable and
/// avoids accidental seed reuse when sweep implementations are rearranged.
pub(crate) fn continuation_seed_for_point(
    policy: ContinuationSeedPolicy,
    configured_seed: &Option<Vec<f64>>,
    previous_accepted_seed: &Option<Vec<f64>>,
) -> Option<Vec<f64>> {
    match policy {
        ContinuationSeedPolicy::PreviousAccepted => previous_accepted_seed.clone(),
        ContinuationSeedPolicy::IndependentPerPoint => configured_seed.clone(),
    }
}

impl EquilibriumLogMoles {
    /// Returns `true` when the current mutable facade has enough prepared
    /// thermochemical data to build the RST symbolic problem directly.
    ///
    /// The canonical typed problem can always carry numeric Gibbs closures, but
    /// the RST adapter still needs prepared searchable thermochemical data for
    /// the symbolic bridge. This helper keeps the default policy honest: RST
    /// Checks whether the solver has the symbolic Gibbs expressions needed by the RST backend.
    ///
    /// The RST backend requires `gibbs_sym` (symbolic `Expr` closures) in addition to
    /// the numeric `GibbsFn` closures. This method returns `true` only when the symbolic
    /// vector is non-empty and its length matches the substance database. Falling back to
    /// the legacy backend when this returns `false` prevents a cryptic RST initialization
    /// failure.
    pub(crate) fn has_rst_symbolic_context(&self) -> bool {
        (self.gibbs_sym.len() == self.subs_data.substances.len() && !self.gibbs_sym.is_empty())
            || !self.subs_data.search_results.is_empty()
            || !self.subs_data.search_states.is_empty()
            || !self.subs_data.therm_map_of_sym.is_empty()
            || !self.subs_data.therm_map_of_fun.is_empty()
    }

    /// Validates a temperature-grid request before it can enter a mutable
    /// coefficient-selection or solve path.
    ///
    /// All three sweep implementations use this gate so an invalid step cannot
    /// Validates temperature range parameters before any sweep setup.
    ///
    /// Checks that:
    /// - All values are finite and strictly positive (T > 0 K).
    /// - `T_end >= T_start` (non-decreasing range).
    /// - `T_step > 0` (strictly positive step).
    /// - The range produces at least one point.
    ///
    /// This prevents division by zero, infinite loops, and nonsensical temperature
    /// inputs before any solver resources are allocated.
    pub(crate) fn validate_temperature_range(
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<(), ReactionExtentError> {
        for (field, value) in [("T_start", T_start), ("T_end", T_end), ("T_step", T_step)] {
            if !value.is_finite() {
                return Err(ReactionExtentError::InvalidProblem {
                    field,
                    message: "temperature-range values must be finite".to_string(),
                });
            }
        }
        if T_step <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "T_step",
                message: "temperature step must be positive".to_string(),
            });
        }
        if T_start >= T_end {
            return Err(ReactionExtentError::InvalidProblem {
                field: "temperature_range",
                message: "T_start must be strictly smaller than T_end".to_string(),
            });
        }
        Ok(())
    }

    /// Builds ordered Gibbs closures from the current `SubsData` coefficients.
    ///
    /// The substance list is the canonical ordering for all solver vectors, so
    /// a missing entry is a typed setup error rather than an `unwrap()` panic.
    pub(crate) fn ordered_gibbs_functions(
        substances: &[String],
        mut functions: HashMap<String, Box<dyn Fn(f64) -> f64 + Send + Sync>>,
    ) -> Result<Vec<Box<dyn Fn(f64) -> f64 + Send + Sync>>, ReactionExtentError> {
        substances
            .iter()
            .map(|substance| {
                functions.remove(substance).ok_or_else(|| {
                    ReactionExtentError::SubsDataError(SubsDataError::MissingData {
                        field: "standard Gibbs closure".to_string(),
                        substance: substance.clone(),
                    })
                })
            })
            .collect()
    }

    /// Builds ordered Gibbs closures from the current `SubsData` coefficients.
    ///
    /// The substance list is the canonical ordering for all solver vectors, so
    /// Builds Gibbs closures for all species in the solver's current substance order.
    ///
    /// Delegates to [`Self::ordered_gibbs_functions`] using the solver's own
    /// `subs_data.substances` list. This is the single-threaded path used by
    /// sequential temperature sweeps and single-point solves.
    pub(crate) fn build_gibbs_functions(&mut self) -> Result<Vec<GibbsFn>, ReactionExtentError> {
        let substances = self.subs_data.substances.clone();
        Self::ordered_gibbs_functions(&substances, self.subs_data.calculate_dG0_fun_one_phase()?)?
            .into_iter()
            .map(|function| Ok(Rc::new(move |temperature: f64| function(temperature)) as GibbsFn))
            .collect()
    }

    /// Builds thread-safe Gibbs closures with the same ordering and missing-data
    /// contract as [`Self::build_gibbs_functions`].
    pub(crate) fn build_parallel_gibbs_functions(
        &mut self,
    ) -> Result<Vec<std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>>, ReactionExtentError> {
        let substances = self.subs_data.substances.clone();
        Self::ordered_gibbs_functions(&substances, self.subs_data.calculate_dG0_fun_one_phase()?)?
            .into_iter()
            .map(|function| {
                Ok(
                    std::sync::Arc::new(move |temperature: f64| function(temperature))
                        as std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>,
                )
            })
            .collect()
    }

    /// Copies the complete accepted-solve evidence into one sweep snapshot.
    pub(crate) fn snapshot_published_solution(
        &self,
        temperature: f64,
    ) -> Result<TemperatureSolveSnapshot, ReactionExtentError> {
        let validation = self.last_validation_report.clone().ok_or_else(|| {
            ReactionExtentError::InvalidProblem {
                field: "temperature_snapshot",
                message: "accepted temperature solve did not publish validation evidence"
                    .to_string(),
            }
        })?;
        let solve_report =
            self.last_solve_report
                .clone()
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "temperature_snapshot",
                    message: "accepted temperature solve did not publish a backend report"
                        .to_string(),
                })?;
        Ok(TemperatureSolveSnapshot {
            temperature,
            log_moles: self.solution.clone(),
            moles: self.moles.clone(),
            validation,
            solve_report,
        })
    }

    /// Builds the legacy mutable engine from one validated canonical problem.
    ///
    /// New callers should construct [`EquilibriumProblem`] first. This adapter
    /// exists while the existing temperature-sweep and phase-control code is
    /// migrated away from its public mutable fields.
    pub fn from_problem(problem: EquilibriumProblem) -> Result<Self, ReactionExtentError> {
        let prepared = PreparedEquilibriumProblem::new(problem)?;
        let (
            species,
            initial_moles,
            initial_log_moles,
            element_composition,
            gibbs,
            phases,
            conditions,
            reaction_basis,
            element_totals,
            species_phase,
        ) = prepared.into_legacy_parts();
        let mut solver = Self::empty();
        solver.subs_data.substances = species;
        solver.n0 = initial_moles;
        solver.initial_guess = Some(initial_log_moles);
        solver.elem_composition = element_composition;
        solver.gibbs = gibbs;
        solver.phases = phases;
        solver.T = conditions.temperature();
        solver.P = conditions.pressure();
        solver.p0 = conditions.reference_pressure();
        solver.stoich_matrix = reaction_basis.reactions.clone();
        solver.reaction_basis = reaction_basis;
        solver.elements_vector = element_totals;
        solver.species_phase = species_phase;
        Ok(solver)
    }

    /// Solves one canonical equilibrium problem with the default solver
    /// configuration.
    ///
    /// This is the explicit construction path the refactor is steering toward:
    /// build a validated problem, optionally adjust typed solver settings, then
    /// perform one solve and return the accepted immutable snapshot.
    pub fn solve_problem(
        problem: EquilibriumProblem,
    ) -> Result<EquilibriumSolution, ReactionExtentError> {
        Self::solve_problem_with(problem, |_| {})
    }

    /// Solves one canonical equilibrium problem after applying a typed solver
    /// configuration hook.
    ///
    /// The hook receives only [`EquilibriumSolverSettings`], so callers can
    /// tune policies and budgets without mutating the already validated
    /// chemistry, phase layout, or thermodynamic conditions.
    pub fn solve_problem_with<F>(
        problem: EquilibriumProblem,
        configure: F,
    ) -> Result<EquilibriumSolution, ReactionExtentError>
    where
        F: FnOnce(&mut EquilibriumSolverSettings),
    {
        let mut solver = Self::from_problem(problem)?;
        configure(&mut solver.solver_settings);
        solver.solve()?;
        solver.accepted_solution()
    }

    /// Returns an immutable snapshot only when the current candidate passed
    /// the backend-independent acceptance gate.
    ///
    /// Legacy mutable fields remain available during migration, but new
    /// callers should use this method rather than combining them manually.
    pub fn accepted_solution(&self) -> Result<EquilibriumSolution, ReactionExtentError> {
        let validation = self.last_validation_report.clone().ok_or_else(|| {
            ReactionExtentError::InvalidCandidate {
                field: "accepted_solution",
                message: "no accepted equilibrium solution is available".to_string(),
            }
        })?;
        let conditions = crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions::new(
            self.T, self.P, self.p0,
        )?;
        EquilibriumSolution::new(
            self.solution.clone(),
            self.moles.clone(),
            conditions,
            validation,
        )
    }

    /// Builds a typed canonical problem snapshot from the current mutable state.
    ///
    /// The snapshot is used by secondary validation paths so they can work on
    /// the same thermodynamic inputs as the main solve without peeking into
    /// solver-only mutable publication fields.
    fn current_problem_snapshot(&self) -> Result<EquilibriumProblem, ReactionExtentError> {
        let initial_log_moles = match &self.initial_guess {
            Some(seed) => LogMolesInitialGuess::new(seed.clone())?,
            None => LogMolesInitialGuess::from_initial_moles(&self.n0)?,
        };
        let conditions = EquilibriumConditions::new(self.T, self.P, self.p0)?;
        EquilibriumProblem::new(
            self.subs_data.substances.clone(),
            self.n0.clone(),
            initial_log_moles,
            self.elem_composition.clone(),
            self.gibbs.clone(),
            self.phases.clone(),
            conditions,
        )
    }

    /// Runs the optional independent K_eq cross-validation boundary.
    ///
    /// This is intentionally separate from the main backend acceptance gate:
    /// the canonical solver still publishes its own acceptance evidence first,
    /// and this secondary pass only decides whether an additional independent
    /// report should be attached or whether a required validation must fail.
    fn run_keq_cross_validation(
        &self,
        candidate: &EquilibriumSolution,
    ) -> Result<Option<EquilibriumConstantCrossValidationStatus>, ReactionExtentError> {
        let mode = self.solver_settings.keq_validation_mode;
        if mode == EquilibriumConstantValidationMode::Off {
            return Ok(None);
        }

        if self.phases.len() != 1 || !matches!(self.phases[0].kind, PhaseActivityModel::IdealGas) {
            let status = EquilibriumConstantCrossValidationStatus::ValidatorNotApplicable {
                message: "independent validation currently supports exactly one ideal-gas phase"
                    .to_string(),
            };
            return match mode {
                EquilibriumConstantValidationMode::Required => {
                    Err(ReactionExtentError::ValidationNotApplicable {
                        path: "equilibrium_constant_cross_validation",
                        message:
                            "independent validation currently supports exactly one ideal-gas phase"
                                .to_string(),
                    })
                }
                EquilibriumConstantValidationMode::Off
                | EquilibriumConstantValidationMode::WhenApplicable => Ok(Some(status)),
            };
        }

        let prepared = PreparedEquilibriumProblem::new(self.current_problem_snapshot()?)?;
        let keq_problem =
            EquilibriumConstantProblem::from_prepared_ideal_gas(&prepared, Default::default())?;
        let solver = EquilibriumConstantSolver {
            mode: match mode {
                EquilibriumConstantValidationMode::Off => EquilibriumConstantSolverMode::Off,
                EquilibriumConstantValidationMode::WhenApplicable => {
                    EquilibriumConstantSolverMode::WhenApplicable
                }
                EquilibriumConstantValidationMode::Required => {
                    EquilibriumConstantSolverMode::Required
                }
            },
            validation_tolerances: EquilibriumConstantValidationTolerances::default(),
            ..Default::default()
        };
        let validator = solver.solve_if_applicable(&keq_problem)?;
        let status = classify_equilibrium_constant_cross_validation(
            &keq_problem,
            Ok(candidate.clone()),
            Ok(validator),
            self.solver_settings.keq_validation_tolerances,
        )?;

        if mode == EquilibriumConstantValidationMode::Required {
            match &status {
                EquilibriumConstantCrossValidationStatus::Compared(report) if report.accepted => {}
                _ => {
                    return Err(ReactionExtentError::InvalidCandidate {
                        field: "equilibrium_constant_cross_validation",
                        message: status.to_string(),
                    });
                }
            }
        }

        Ok(Some(status))
    }

    /// Creates a solver with default settings and an empty substance database.
    ///
    /// The solver starts with `T = 273.15 K`, `P = 101325 Pa`, `p0 = 101325 Pa`,
    /// no species, no elements, and default solver parameters (LM backend, no scaling).
    /// Call [`set_problem()`](Self::set_problem) or populate fields directly before solving.
    pub fn new() -> Self {
        let ReactionBasis0 = ReactionBasis {
            rank: 0,
            num_reactions: 0,
            reactions: DMatrix::zeros(0, 0),
        };
        Self {
            subs_data: SubsData::new(),
            elem_composition: DMatrix::zeros(0, 0),
            reaction_basis: ReactionBasis0,
            initial_guess: None,
            P: 101325.0,
            T: 273.15,
            stoich_matrix: DMatrix::zeros(0, 0),
            n0: Vec::new(),
            gibbs: Vec::new(),
            gibbs_sym: Vec::new(),
            phases: Vec::new(),
            species_phase: Vec::new(),
            phase_active_mask: Vec::new(),
            solution: Vec::new(),
            last_validation_report: None,
            last_solve_report: None,
            last_keq_validation_status: None,
            last_phase_control_report: None,
            list_of_failed_T: Vec::new(),
            temperature_failures: Vec::new(),
            temperature_solutions: Vec::new(),
            map_of_moles_for_each_substance: HashMap::new(),
            moles_for_T_range: Vec::new(),
            moles: Vec::new(),
            elements_vector: Vec::new(),
            loglevel: None,
            species_eps: DEFAULT_TRACE_MOLE_FLOOR,
            substate_eps: DEFAULT_TRACE_MOLE_FLOOR,
            solver_settings: EquilibriumSolverSettings::default(),
            phase_manager: PhaseManager::default(),
            p0: 101325.0,
        }
    }

    /// Creates a solver with default settings and a completely empty database.
    ///
    /// Identical to [`new()`](Self::new) but uses `SubsData::empty()` instead of
    /// `SubsData::new()`. Both constructors produce the same initial state for
    /// practical purposes. Provided for semantic clarity when the caller explicitly
    /// wants an empty database.
    pub fn empty() -> Self {
        let ReactionBasis0 = ReactionBasis {
            rank: 0,
            num_reactions: 0,
            reactions: DMatrix::zeros(0, 0),
        };
        Self {
            subs_data: SubsData::empty(),
            elem_composition: DMatrix::zeros(0, 0),
            reaction_basis: ReactionBasis0,
            initial_guess: None,
            P: 101325.0,
            T: 273.15,
            stoich_matrix: DMatrix::zeros(0, 0),
            n0: Vec::new(),
            gibbs: Vec::new(),
            gibbs_sym: Vec::new(),
            phases: Vec::new(),
            species_phase: Vec::new(),
            phase_active_mask: Vec::new(),
            solution: Vec::new(),
            last_validation_report: None,
            last_solve_report: None,
            last_keq_validation_status: None,
            last_phase_control_report: None,
            list_of_failed_T: Vec::new(),
            temperature_failures: Vec::new(),
            temperature_solutions: Vec::new(),
            map_of_moles_for_each_substance: HashMap::new(),
            moles_for_T_range: Vec::new(),
            moles: Vec::new(),
            elements_vector: Vec::new(),
            loglevel: None,
            species_eps: DEFAULT_TRACE_MOLE_FLOOR,
            substate_eps: DEFAULT_TRACE_MOLE_FLOOR,
            solver_settings: EquilibriumSolverSettings::default(),
            phase_manager: PhaseManager::default(),
            p0: 101325.0,
        }
    }
    /// Computes the stoichiometric matrix `ν` from the SVD reaction basis.
    ///
    /// The stoichiometric matrix has shape `(species × reactions)` where each column
    /// is one independent reaction direction from the SVD nullspace of the element
    /// composition matrix. This must be called after `elem_composition` is populated
    /// and before any solver backend is invoked.
    pub fn create_stoich_matrix(&mut self) -> Result<(), ReactionExtentError> {
        self.solver_settings.validate()?;
        let reaction_basis = compute_reaction_basis(
            &self.elem_composition,
            self.solver_settings.solver_params.tol,
        )?;
        let stoich_matrix = reaction_basis.reactions.clone();
        let elements_vector = compute_element_totals(&self.elem_composition, &self.n0)?;
        let species_phase = species_to_phase_map(&self.phases, self.n0.len())?;

        self.reaction_basis = reaction_basis;
        self.stoich_matrix = stoich_matrix;
        self.elements_vector = elements_vector.data.as_vec().clone();
        self.species_phase = species_phase;
        Ok(())
    }

    /// Stores the caller's preferred log level for this solver request.
    ///
    /// A numerical solver must not initialize a process-global logger: GUI,
    /// library, and test hosts own that responsibility. The setting is kept as
    /// request metadata while legacy workflow constructors migrate to typed
    /// Sets the log level for diagnostic output and returns the solver for chaining.
    ///
    /// Pass `None` to disable logging, or `Some("info")` / `Some("debug")` for
    /// progressively more verbose output. This is a builder-pattern convenience
    /// wrapper around the mutable `loglevel` field.
    pub fn with_loglevel(&mut self, loglevel: Option<&str>) -> &mut Self {
        self.loglevel = loglevel.map(str::to_string);
        self
    }
    /// Sets a legacy log-mole initial iterate.
    ///
    /// The values are `ln(n_i)`, not physical mole numbers. New code should
    /// prefer `LogMolesInitialGuess` through `EquilibriumProblem`.
    /// Replaces the log-mole initial iterate after validating it immediately.
    ///
    /// This legacy mutable entry point follows the same finite-value and
    /// dimension contract as [`LogMolesInitialGuess`], rather than delaying a
    /// Overrides the default initial guess with an explicit log-mole seed.
    ///
    /// The seed is validated (all values must be finite) before it replaces the
    /// current `initial_guess`. This allows callers to provide a warm-start from
    /// a previous solve or an external estimate, avoiding the default trace-floor
    /// initialization.
    pub fn set_initial_guess(&mut self, guess: Vec<f64>) -> Result<(), ReactionExtentError> {
        let guess = LogMolesInitialGuess::new(guess)?;
        let expected = self.subs_data.substances.len();
        if expected > 0 && guess.as_slice().len() != expected {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial log-mole guess has {} entries, expected {expected}",
                guess.as_slice().len(),
            )));
        }
        self.initial_guess = Some(guess.into_inner());
        Ok(())
    }

    /// Replaces the default backend policy with an explicit [`SolverPolicy`].
    ///
    /// When set, this policy determines the backend order (single or cascade)
    /// instead of the legacy `solver` field. Pass a `Cascade` with specific
    /// backends to control fallback behavior.
    pub fn set_solver_policy(&mut self, policy: SolverPolicy) {
        self.solver_settings.solver_policy = Some(policy);
    }

    /// Clears the explicit solver policy, reverting to the legacy default order.
    ///
    /// After calling this, the backend order is determined by the `solver` field
    /// and the historical fallback chain (preferred → LM → NR → TR).
    pub fn clear_solver_policy(&mut self) {
        self.solver_settings.solver_policy = None;
    }

    /// Sets explicit iteration and attempt limits for the solver cascade.
    ///
    /// The budget controls how many backends may start (`max_attempts`), how many
    /// iterations each backend receives (`max_iterations_per_attempt`), and the
    /// total iteration budget across all backends (`max_total_iterations`).
    pub fn set_solver_budget(&mut self, budget: SolverCascadeBudget) {
        self.solver_settings.solver_budget = Some(budget);
    }

    /// Clears the explicit cascade budget, reverting to the policy-derived default.
    ///
    /// After calling this, the budget is computed from `SolverParams.max_iter`
    /// and the number of backends in the current policy.
    pub fn clear_solver_budget(&mut self) {
        self.solver_settings.solver_budget = None;
    }

    /// Selects whether a temperature sweep continues from accepted points or
    /// Sets the warm-start behavior for sequential temperature sweeps.
    ///
    /// - `PreviousAccepted` (default): each temperature point starts from the
    ///   previous point's solution (faster convergence for smooth grids).
    /// - `IndependentPerPoint`: each point starts from the configured trace seed
    ///   (required for fully parallel sweeps).
    pub fn set_continuation_seed_policy(&mut self, policy: ContinuationSeedPolicy) {
        self.solver_settings.continuation_seed_policy = policy;
    }

    /// Clears every published solve artifact after a new mutable problem state
    /// has been validated and applied.
    pub(crate) fn clear_published_state(&mut self) {
        self.solution.clear();
        self.moles.clear();
        self.last_validation_report = None;
        self.last_solve_report = None;
        self.last_keq_validation_status = None;
        self.last_phase_control_report = None;
        self.list_of_failed_T.clear();
        self.temperature_failures.clear();
        self.temperature_solutions.clear();
        self.moles_for_T_range.clear();
        self.map_of_moles_for_each_substance.clear();
    }

    /// Validates the legacy mutable problem shape before the state is staged.
    pub(crate) fn validate_mutable_problem_shape(&self, n0: &[f64]) -> Result<(), ReactionExtentError> {
        if n0.iter().any(|moles| !moles.is_finite() || *moles < 0.0) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_moles",
                message: "initial moles must be finite and non-negative".to_string(),
            });
        }
        if !self.subs_data.substances.is_empty() && n0.len() != self.subs_data.substances.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial moles have {} entries but SubsData contains {} substances",
                n0.len(),
                self.subs_data.substances.len(),
            )));
        }
        if let Some(initial_guess) = self.initial_guess.as_ref() {
            if initial_guess.len() != n0.len() {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "initial log-mole guess has {} entries, but the new problem has {} species",
                    initial_guess.len(),
                    n0.len(),
                )));
            }
        }
        Ok(())
    }

    /// Stages the mutable equilibrium-system inputs in local values first.
    pub(crate) fn stage_equilibrium_system(
        &self,
        n0: &[f64],
    ) -> Result<
        (
            SubsData,
            DMatrix<f64>,
            ReactionBasis,
            DMatrix<f64>,
            Vec<f64>,
            Vec<usize>,
        ),
        ReactionExtentError,
    > {
        self.solver_settings.validate()?;
        let mut staged_subs_data = self.subs_data.clone();
        staged_subs_data
            .search_substances()
            .map_err(ReactionExtentError::SubsDataError)?;
        staged_subs_data
            .parse_all_thermal_coeffs()
            .map_err(ReactionExtentError::SubsDataError)?;
        staged_subs_data
            .calculate_elem_composition_and_molar_mass(None)
            .map_err(ReactionExtentError::SubsDataError)?;
        let element_composition = staged_subs_data
            .element_composition_matrix()
            .cloned()
            .ok_or_else(|| {
                ReactionExtentError::SubsDataError(SubsDataError::MissingData {
                    field: "element composition matrix after calculation".to_string(),
                    substance: "all_substances".to_string(),
                })
            })?;

        if n0.len() != element_composition.nrows() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial moles have {} entries but the staged element matrix has {} species",
                n0.len(),
                element_composition.nrows(),
            )));
        }

        let reaction_basis =
            compute_reaction_basis(&element_composition, self.solver_settings.solver_params.tol)?;
        let elements_vector = compute_element_totals(&element_composition, n0)?;
        let species_phase = species_to_phase_map(&self.phases, n0.len())?;
        let stoich_matrix = reaction_basis.reactions.clone();

        Ok((
            staged_subs_data,
            element_composition,
            reaction_basis,
            stoich_matrix,
            elements_vector.data.as_vec().clone(),
            species_phase,
        ))
    }

    /// Resolves the active log-mole seed for the current mutable problem.
    ///
    /// If the caller did not set an explicit seed, the helper derives one from
    /// the physical initial moles using the canonical trace floor. That keeps
    /// the fallback in one place and prevents the same seed contract from
    /// being re-encoded in `solve`, sweep setup, and tests.
    pub(crate) fn resolved_initial_guess(&self) -> Result<(Vec<f64>, bool), ReactionExtentError> {
        match &self.initial_guess {
            Some(seed) => {
                let guess = LogMolesInitialGuess::new(seed.clone())?;
                let expected = self.subs_data.substances.len();
                if expected > 0 && guess.as_slice().len() != expected {
                    return Err(ReactionExtentError::DimensionMismatch(format!(
                        "initial log-mole guess has {} entries, expected {expected}",
                        guess.as_slice().len(),
                    )));
                }
                Ok((guess.into_inner(), false))
            }
            None => Ok((
                LogMolesInitialGuess::from_moles_with_policy(
                    &self.n0,
                    self.solver_settings.trace_seed_policy,
                )?
                .into_inner(),
                true,
            )),
        }
    }

    /// Builds the nonlinear residual/Jacobian contract used by the solver
    /// cascade.
    ///
    /// Keeping this in one place ensures the scalar solve and all sweep
    /// variants evaluate the same equations in the same order.
    #[allow(dead_code)]
    pub(crate) fn build_solve_contract(&mut self) -> Result<EquilibriumSolveContract, ReactionExtentError> {
        let subs_eps = self.substate_eps;
        let phase_eps = self.species_eps;
        let stoich = Rc::new(self.stoich_matrix.clone());
        let species_phase = Rc::new(self.species_phase.clone());
        let delta_n = Rc::new(reaction_phase_stoichiometry(
            &self.stoich_matrix,
            &self.phases,
        ));
        let elem_composition = self.elem_composition.clone();
        let elements_vector = self.elements_vector.clone();

        let f_raw: Rc<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>> =
            Rc::from(equilibrium_logmole_residual(
                (*stoich).clone(),
                elem_composition.clone(),
                elements_vector.clone(),
                self.gibbs.clone(),
                self.phases.clone(),
                self.T,
                self.P,
                self.p0,
                (*species_phase).clone(),
                subs_eps,
                phase_eps,
            )?);

        let scale = equilibrium_scaling(
            &stoich,
            &self.elem_composition,
            &self.gibbs,
            &elements_vector,
            self.T,
        )?;
        let f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>> =
            if self.solver_settings.scaling_flag {
                let f_raw = f_raw.clone();
                let residual_scale = scale.clone();
                Box::new(move |values| scale_residual_rows(f_raw(values)?, &residual_scale))
            } else {
                let f_raw = f_raw.clone();
                Box::new(move |values| f_raw(values))
            };

        let ordered_backends = self
            .solver_settings
            .solver_policy
            .clone()
            .unwrap_or_else(|| {
                if self.has_rst_symbolic_context() {
                    SolverPolicy::rusted_scithe_default()
                } else {
                    SolverPolicy::legacy_default(self.solver_settings.solver)
                }
            })
            .ordered_backends();
        let uses_legacy_backend = ordered_backends
            .iter()
            .any(|backend| matches!(backend, SolverBackend::Legacy(_)));
        let uses_rst_backend = ordered_backends
            .iter()
            .any(|backend| matches!(backend, SolverBackend::RustedSciThe(_)));

        let phases_len = self.phases.len();
        let scaling_flag = self.solver_settings.scaling_flag;
        let j = if uses_legacy_backend {
            let stoich = stoich.clone();
            let species_phase = species_phase.clone();
            let delta_n = delta_n.clone();
            let elem_composition = elem_composition.clone();
            let j_raw = Box::new(move |y: &[f64]| {
                equilibrium_logmole_jacobian(
                    y,
                    &stoich,
                    &elem_composition,
                    &species_phase,
                    delta_n.as_ref(),
                    phases_len,
                    subs_eps,
                    phase_eps,
                )
            })
                as Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>;
            Some(if scaling_flag {
                scaled_jacobian(j_raw, scale)
            } else {
                j_raw
            })
        } else {
            None
        };

        let rst_problem = if uses_rst_backend {
            Some(prepare_rst_symbolic_problem(self)?)
        } else {
            None
        };

        Ok(EquilibriumSolveContract {
            f_raw,
            f,
            j,
            rst_problem,
        })
    }

    /// Captures the mutable sweep-relevant state without borrowing the main
    /// solver object in a rayon closure.
    pub(crate) fn temperature_worker_seed(&self) -> TemperatureWorkerSeed {
        TemperatureWorkerSeed {
            subs_data: self.subs_data.clone(),
            elem_composition: self.elem_composition.clone(),
            reaction_basis: self.reaction_basis.clone(),
            stoich_matrix: self.stoich_matrix.clone(),
            n0: self.n0.clone(),
            phases: self.phases.clone(),
            species_phase: self.species_phase.clone(),
            elements_vector: self.elements_vector.clone(),
            loglevel: self.loglevel.clone(),
            species_eps: self.species_eps,
            substate_eps: self.substate_eps,
            solver_settings: self.solver_settings.clone(),
            P: self.P,
            p0: self.p0,
        }
    }

    /// Finds contiguous temperature intervals where the coefficient cache can
    /// be reused without recalculating Gibbs polynomials.
    ///
    /// Both sequential and parallel sweeps use the same interval discovery
    /// logic so a coefficient update is handled once, in one place.
    pub(crate) fn build_temperature_ranges(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, f64)>, ReactionExtentError> {
        let mut temp_ranges = Vec::new();
        let mut T = T_start;

        while T < T_end {
            let user_subs = &mut self.subs_data;
            let vec_of_coeffs = user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(T),
                    })
                })?;

            if vec_of_coeffs.len() != 0 || self.gibbs.len() == 0 {
                let range_start = T;
                let mut T_next = T + T_step;
                while T_next < T_end {
                    let changed = self
                        .subs_data
                        .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T_next)
                        .map_err(|_| {
                            ReactionExtentError::SubsDataError(
                                SubsDataError::CoefficientExtractionFailed {
                                    substance: "all_subs".to_string(),
                                    temperature: Some(T_next),
                                },
                            )
                        })?;
                    if !changed.is_empty() {
                        break;
                    }
                    T_next += T_step;
                }
                temp_ranges.push((range_start, T_next));
                T = T_next;
            } else {
                T += T_step;
            }
        }

        Ok(temp_ranges)
    }

    /// Builds the precomputed Gibbs-function cache for each temperature band.
    ///
    /// This is intentionally separate from interval discovery so the code can
    /// reuse the same range plan for sequential and parallel execution.
    pub(crate) fn build_temperature_gibbs_cache(
        &mut self,
        temp_ranges: &[(f64, f64)],
    ) -> Result<
        std::collections::HashMap<
            usize,
            (
                Vec<std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>>,
                f64,
                f64,
            ),
        >,
        ReactionExtentError,
    > {
        let mut gibbs_cache = std::collections::HashMap::new();

        for (range_id, (T_start_range, T_end_range)) in temp_ranges.iter().enumerate() {
            self.subs_data
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(*T_start_range)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(*T_start_range),
                    })
                })?;

            let g_vec = self.build_parallel_gibbs_functions()?;
            gibbs_cache.insert(range_id, (g_vec, *T_start_range, *T_end_range));
        }

        Ok(gibbs_cache)
    }

    /// Maps each temperature step to the Gibbs-cache band that covers it.
    pub(crate) fn build_temperature_point_index(
        temp_ranges: &[(f64, f64)],
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Vec<(f64, usize)> {
        let mut temp_range_pairs = Vec::new();
        let mut T = T_start;

        while T < T_end {
            if let Some((range_id, _)) = temp_ranges
                .iter()
                .enumerate()
                .find(|(_, (t_start, t_end))| T >= *t_start && T < *t_end)
            {
                temp_range_pairs.push((T, range_id));
            }
            T += T_step;
        }

        temp_range_pairs
    }

    /// Solves one temperature point inside a sweep using a captured canonical
    /// seed and one precomputed Gibbs-function set.
    ///
    /// Both parallel sweep variants call this helper so the worker setup,
    /// backend contract, and snapshot publication all stay identical.
    pub(crate) fn solve_temperature_point_from_seed(
        worker_seed: &TemperatureWorkerSeed,
        temperature: f64,
        gibbs_vec: Vec<std::sync::Arc<dyn Fn(f64) -> f64 + Send + Sync>>,
        initial_guess: Option<Vec<f64>>,
        point_index: f64,
    ) -> Result<(TemperatureSolveSnapshot, f64), (f64, ReactionExtentError)> {
        let now = Instant::now();
        let mut local_solver = EquilibriumLogMoles::empty();

        let gibbs_as_rc: Vec<GibbsFn> = gibbs_vec
            .iter()
            .map(|arc_fn| {
                let fn_clone = arc_fn.clone();
                Rc::new(move |T: f64| fn_clone(T)) as GibbsFn
            })
            .collect();

        worker_seed.apply(&mut local_solver, temperature, gibbs_as_rc, initial_guess);

        match local_solver.solve() {
            Ok(()) => local_solver
                .snapshot_published_solution(temperature)
                .map(|snapshot| (snapshot, now.elapsed().as_millis() as f64))
                .map_err(|error| (point_index, error)),
            Err(e) => Err((point_index, e)),
        }
    }

    /// Replaces the mutable problem state after validating it as one unit.
    ///
    /// Previously published solve artifacts are cleared only after validation
    /// succeeds, so a rejected update cannot leave a half-reconfigured solver.
    pub fn set_problem(
        &mut self,
        n0: Vec<f64>,
        phases: Vec<Phase>,
        P: f64,
    ) -> Result<(), ReactionExtentError> {
        if !P.is_finite() || P <= 0.0 {
            return Err(ReactionExtentError::InvalidConditions {
                parameter: "pressure",
                value: P,
            });
        }
        self.validate_mutable_problem_shape(&n0)?;
        let species_phase = species_to_phase_map(&phases, n0.len())?;

        self.n0 = n0;
        self.phases = phases;
        self.P = P;
        self.species_phase = species_phase;
        self.phase_active_mask.clear();
        self.clear_published_state();
        Ok(())
    }
    /// Applies an explicit initial-mole map without accepting misspelled names
    /// or publishing a partially updated vector.
    pub fn set_n0_from_non_zero_map(
        &mut self,
        map_of_nonzero_moles: HashMap<String, f64>,
    ) -> Result<(), ReactionExtentError> {
        let subs = self.subs_data.substances.clone();
        for (substance, moles) in &map_of_nonzero_moles {
            if !subs.iter().any(|known| known == substance) {
                return Err(ReactionExtentError::SubsDataError(
                    SubsDataError::SubstanceNotFound(substance.clone()),
                ));
            }
            if !moles.is_finite() || *moles < 0.0 {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_moles",
                    message: format!(
                        "initial amount for '{substance}' must be finite and non-negative"
                    ),
                });
            }
        }
        let mut n0 = vec![1e-8; subs.len()];
        for (i, subs_i) in subs.iter().enumerate() {
            if let Some(moles_i) = map_of_nonzero_moles.get(subs_i) {
                n0[i] = *moles_i;
            }
        }
        self.n0 = n0;
        self.phase_active_mask.clear();
        self.clear_published_state();
        Ok(())
    }
    /// Creates a fully prepared equilibrium system as one transactional update.
    ///
    /// Substance lookup, coefficient parsing, elemental composition, and the
    /// reaction basis are built in local values first. The solver publishes the
    /// new state only when every stage succeeds, avoiding a mixture of old and
    /// new matrices after an incomplete setup attempt.
    pub fn create_equilibrium_system(&mut self) -> Result<(), ReactionExtentError> {
        let (
            staged_subs_data,
            element_composition,
            reaction_basis,
            stoich_matrix,
            elements_vector,
            species_phase,
        ) = self.stage_equilibrium_system(&self.n0)?;
        self.subs_data = staged_subs_data;
        self.elem_composition = element_composition;
        self.reaction_basis = reaction_basis;
        self.stoich_matrix = stoich_matrix;
        self.elements_vector = elements_vector;
        self.species_phase = species_phase;
        self.phase_active_mask.clear();

        Ok(())
    }
    pub fn create_symbolic_system(&mut self) -> Result<(), ReactionExtentError> {
        // creating variables
        let stoich = self.stoich_matrix.clone();
        let vars = Expr::IndexedVars(stoich.nrows(), "y").0;
        let vars_Str: Vec<String> = vars.iter().map(|y| y.to_string().clone()).collect();
        //let vars_str: Vec<&str> = vars_Str.iter().map(|yi| yi.as_str()).collect();
        let mut jac = Jacobian::new();
        // creating vector of symbolic dG0
        let user_subs = &mut self.subs_data;
        user_subs
            .calculate_therm_map_of_sym()
            .map_err(|e| ReactionExtentError::SubsDataError(e))?;
        let map_of_functions = user_subs.calculate_dG0_sym_one_phase()?;
        let mut g_vec: Vec<Expr> = Vec::new();
        for subs in user_subs.substances.clone() {
            let dg0 = map_of_functions.get(&subs).ok_or_else(|| {
                ReactionExtentError::SubsDataError(SubsDataError::MissingData {
                    field: "standard Gibbs symbolic expression".to_string(),
                    substance: subs.clone(),
                })
            })?;
            g_vec.push(dg0.clone());
        }
        // vector of elements
        // let elements_vector = compute_element_totals(&self.elem_composition, &self.n0);
        // self.elements_vector = elements_vector.data.as_vec().clone();
        let elements_vector = self.elements_vector.clone();
        // creating vector of symbolic equations for equilibrium
        let vec_of_functions = multiphase_equilibrium_residual_generator_sym(
            self.stoich_matrix.clone(),
            self.elem_composition.clone(),
            elements_vector.clone(),
            g_vec,
            self.phases.clone(),
            self.P,
            self.p0,
        )?;
        jac.set_vector_of_functions(vec_of_functions);
        jac.vector_of_variables = vars;
        jac.variable_string = vars_Str;
        jac.parameters_string = vec!["T".to_string()];
        jac.calc_jacobian();

        Ok(())
    }
    pub fn solve_for_T_range(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        Self::validate_temperature_range(T_start, T_end, T_step)?;
        let now = Instant::now();
        let mut T = T_start;
        let mut moles_for_T_range = Vec::new();
        let mut failures = Vec::new();
        let mut snapshots = Vec::new();
        let configured_seed = self.initial_guess.clone();

        while T < T_end {
            info!("\n solving equilibrium for T = {} \n", T);
            self.initial_guess = continuation_seed_for_point(
                self.solver_settings.continuation_seed_policy,
                &configured_seed,
                &self.initial_guess,
            );
            // Build Gibbs functions while holding a short-lived mutable borrow to subs_data
            let user_subs = &mut self.subs_data;
            let vec_of_coeffs = user_subs
                .extract_coeffs_if_current_coeffs_not_valid_for_all_subs(T)
                .map_err(|_| {
                    ReactionExtentError::SubsDataError(SubsDataError::CoefficientExtractionFailed {
                        substance: "all_subs".to_string(),
                        temperature: Some(T),
                    })
                })?;
            if vec_of_coeffs.len() != 0 || self.gibbs.len() == 0 {
                info!("\n coefficients outdated ... \n renewing pipeline ...");
                self.gibbs = self.build_gibbs_functions()?;
            }
            self.T = T; // Temperature in K
            match self.solve() {
                Ok(()) => {}
                Err(e) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    failures.push(temperature_failure(T, &e));
                    T += T_step;
                    continue;
                }
            }

            let snapshot = self.snapshot_published_solution(T)?;
            if self.solver_settings.continuation_seed_policy
                == ContinuationSeedPolicy::PreviousAccepted
            {
                self.initial_guess = Some(snapshot.log_moles.clone());
            }
            moles_for_T_range.push((T, snapshot.moles.clone()));
            snapshots.push(snapshot);
            T += T_step;
        }
        let mole_table = build_mole_series_table(&self.subs_data.substances, &moles_for_T_range)?;
        self.moles_for_T_range = moles_for_T_range.clone();
        self.list_of_failed_T = failures.iter().map(|failure| failure.temperature).collect();
        self.temperature_failures = failures;
        self.temperature_solutions = snapshots;
        self.map_of_moles_for_each_substance = mole_table;
        info!(
            "total elapsed time non-parallel solver{} ms",
            now.elapsed().as_millis()
        );
        Ok(Vec::new())
    }

    /// parallel solver
    pub fn solve_for_T_range_par(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        Self::validate_temperature_range(T_start, T_end, T_step)?;
        use rayon::prelude::*;
        use std::sync::Arc;
        let now = Instant::now();

        let temp_ranges = self.build_temperature_ranges(T_start, T_end, T_step)?;
        let step2 = now.elapsed().as_millis();
        let gibbs_cache = self.build_temperature_gibbs_cache(&temp_ranges)?;
        let step3 = now.elapsed().as_millis();
        let temp_range_pairs =
            Self::build_temperature_point_index(&temp_ranges, T_start, T_end, T_step);
        let step4 = now.elapsed().as_millis();
        let continuation_seed_policy = self.solver_settings.continuation_seed_policy;
        let configured_seed = self.initial_guess.clone();
        let gibbs_cache = Arc::new(gibbs_cache);
        let worker_seed = self.temperature_worker_seed();
        info!("time for closure creation {:?}", now.elapsed());
        info!("closure stages: {step2}, {step3}, {step4}");
        let now = Instant::now();
        // Chunk the temperature tasks to reduce per-task scheduling overhead.
        // Compute chunk size based on Rayon thread count so work is balanced.
        let num_threads = rayon::current_num_threads();
        let chunk_size = std::cmp::max(1, (temp_range_pairs.len() + num_threads - 1) / num_threads);

        let chunks: Vec<Vec<(f64, usize)>> = temp_range_pairs
            .chunks(chunk_size)
            .map(|c| c.to_vec())
            .collect();

        let mut chunk_results: Vec<
            Vec<Result<(TemperatureSolveSnapshot, f64), (f64, ReactionExtentError)>>,
        > = Vec::new();

        for chunk in chunks {
            let initial_guess = continuation_seed_for_point(
                continuation_seed_policy,
                &configured_seed,
                &self.initial_guess,
            );
            info!("starting temperature chunk with seed {:?}", initial_guess);
            let local_results: Vec<
                Result<(TemperatureSolveSnapshot, f64), (f64, ReactionExtentError)>,
            > = chunk
                .par_iter()
                .map(|(T, range_id)| {
                    // Retrieve the Gibbs functions for this range
                    let gibbs_vec = if let Some((g_vec, _, _)) = gibbs_cache.get(range_id) {
                        g_vec.clone()
                    } else {
                        return Err((
                            *T,
                            ReactionExtentError::InvalidProblem {
                                field: "gibbs_cache",
                                message: format!(
                                    "no precomputed Gibbs functions for range {range_id}"
                                ),
                            },
                        ));
                    };
                    Self::solve_temperature_point_from_seed(
                        &worker_seed,
                        *T,
                        gibbs_vec,
                        initial_guess.clone(),
                        *T,
                    )
                })
                .collect();
            let last_solution = local_results.iter().rev().find_map(|result| {
                result
                    .as_ref()
                    .ok()
                    .map(|(snapshot, _)| snapshot.log_moles.clone())
            });
            let time = local_results
                .iter()
                .map(|res| match res {
                    Ok((_, t)) => *t,
                    Err(_) => 0.0,
                })
                .sum::<f64>();
            info!("time for temperature chunk {:?} ms", time);
            if continuation_seed_policy == ContinuationSeedPolicy::PreviousAccepted {
                if let Some(last_solution) = last_solution {
                    // Use only an accepted local result as the next chunk's seed.
                    self.initial_guess = Some(last_solution);
                } else {
                    warn!(
                        "temperature chunk produced no accepted solution; preserving the previous seed"
                    );
                }
            }
            chunk_results.push(local_results);
        }

        // Flatten chunk results into a single results Vec
        let results: Vec<Result<(TemperatureSolveSnapshot, f64), (f64, ReactionExtentError)>> =
            chunk_results.into_iter().flatten().collect();

        // Step 4: Process results
        let mut moles_for_T_range = Vec::new();
        let mut failed_temps = Vec::new();
        let mut failures = Vec::new();
        let mut snapshots = Vec::new();

        for result in results {
            match result {
                Ok((snapshot, _time)) => {
                    moles_for_T_range.push((snapshot.temperature, snapshot.moles.clone()));
                    snapshots.push(snapshot);
                }
                Err((T, e)) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    failed_temps.push(T);
                    failures.push(temperature_failure(T, &e));
                }
            }
        }

        moles_for_T_range.sort_by(|a, b| a.0.total_cmp(&b.0));
        snapshots.sort_by(|a, b| a.temperature.total_cmp(&b.temperature));

        let mole_table = build_mole_series_table(&self.subs_data.substances, &moles_for_T_range)?;
        self.moles_for_T_range = moles_for_T_range.clone();
        self.list_of_failed_T = failed_temps;
        self.temperature_failures = failures;
        self.temperature_solutions = snapshots;
        self.map_of_moles_for_each_substance = mole_table;
        info!("time for parallel solving {:?}", now.elapsed());
        Ok(moles_for_T_range)
    }

    /// Solves every temperature point independently in parallel.
    ///
    /// This mode deliberately does not apply [`ContinuationSeedPolicy`]: a
    /// previous accepted point is unavailable while tasks run concurrently.
    pub fn solve_for_T_range_par2(
        &mut self,
        T_start: f64,
        T_end: f64,
        T_step: f64,
    ) -> Result<Vec<(f64, Vec<f64>)>, ReactionExtentError> {
        Self::validate_temperature_range(T_start, T_end, T_step)?;
        use rayon::prelude::*;
        use std::sync::Arc;
        let now = Instant::now();
        let temp_ranges = self.build_temperature_ranges(T_start, T_end, T_step)?;
        let step2 = now.elapsed().as_millis();
        let gibbs_cache = self.build_temperature_gibbs_cache(&temp_ranges)?;
        let step3 = now.elapsed().as_millis();
        let temp_range_pairs =
            Self::build_temperature_point_index(&temp_ranges, T_start, T_end, T_step);
        let step4 = now.elapsed().as_millis();
        let gibbs_cache = Arc::new(gibbs_cache);
        let worker_seed = self.temperature_worker_seed();
        info!("time for closure creation {:?}", now.elapsed());
        info!("closure stages: {step2}, {step3}, {step4}");
        let results: Vec<Result<TemperatureSolveSnapshot, (f64, ReactionExtentError)>> =
            temp_range_pairs
                .par_iter()
                .map(|(T, range_id)| {
                    let gibbs_vec = if let Some((g_vec, _, _)) = gibbs_cache.get(range_id) {
                        g_vec.clone()
                    } else {
                        return Err((
                            *T,
                            ReactionExtentError::InvalidProblem {
                                field: "gibbs_cache",
                                message: format!(
                                    "no precomputed Gibbs functions for range {range_id}"
                                ),
                            },
                        ));
                    };

                    Self::solve_temperature_point_from_seed(&worker_seed, *T, gibbs_vec, None, *T)
                        .map(|(snapshot, _elapsed_ms)| snapshot)
                })
                .collect();

        // Step 4: Process results
        let mut moles_for_T_range = Vec::new();
        let mut failed_temps = Vec::new();
        let mut failures = Vec::new();
        let mut snapshots = Vec::new();

        for result in results {
            match result {
                Ok(snapshot) => {
                    moles_for_T_range.push((snapshot.temperature, snapshot.moles.clone()));
                    snapshots.push(snapshot);
                }
                Err((T, e)) => {
                    error!("Equilibrium calculation failed at T = {} K: {:?}", T, e);
                    failed_temps.push(T);
                    failures.push(temperature_failure(T, &e));
                }
            }
        }

        moles_for_T_range.sort_by(|a, b| a.0.total_cmp(&b.0));
        snapshots.sort_by(|a, b| a.temperature.total_cmp(&b.temperature));

        let mole_table = build_mole_series_table(&self.subs_data.substances, &moles_for_T_range)?;
        self.moles_for_T_range = moles_for_T_range.clone();
        self.list_of_failed_T = failed_temps;
        self.temperature_failures = failures;
        self.temperature_solutions = snapshots;
        self.map_of_moles_for_each_substance = mole_table;
        info!("time for parallel solving {:?}", now.elapsed());
        Ok(moles_for_T_range)
    }

    /// Publishes a column-oriented sweep table only after every result row has
    /// been checked against the canonical substance ordering.
    pub fn map_of_moles_for_each_substance(&mut self) -> Result<(), ReactionExtentError> {
        let published =
            build_mole_series_table(&self.subs_data.substances, &self.moles_for_T_range)?;
        self.map_of_moles_for_each_substance = published;
        Ok(())
    }

    /// Builds a typed postprocessing view from the latest solved temperature sweep.
    ///
    /// The solver keeps the solved points as the source of truth. This helper
    /// only packages them for plotting, export, or GUI resampling.
    pub fn temperature_postprocessing(
        &self,
        policy: &TemperaturePostprocessingPolicy,
    ) -> Result<TemperaturePostprocessingResult, ReactionExtentError> {
        postprocess_temperature_series(
            self.subs_data.substances.clone(),
            &self.moles_for_T_range,
            policy,
        )
    }

    /// Builds a display table for the currently published sweep results.
    ///
    /// The equilibrium engine deliberately does not print it: callers choose
    /// whether to render the table in a terminal, GUI, file, or test report.
    pub fn moles_table(&self) -> Table {
        let mut table = Table::new();

        // Get the keys (headers) and sort them for consistent order
        let mut keys: Vec<&String> = self.map_of_moles_for_each_substance.keys().collect();
        keys.sort();

        // Create header row
        let header_row: Vec<Cell> = keys.iter().map(|k| Cell::new(k)).collect();
        table.add_row(Row::new(header_row));

        // Determine the number of rows (length of the vectors)
        let num_rows = if let Some(first_vec) = self.map_of_moles_for_each_substance.values().next()
        {
            first_vec.len()
        } else {
            0
        };

        // Add data rows
        for row_idx in 0..num_rows {
            let mut row: Vec<Cell> = Vec::new();
            for key in &keys {
                if let Some(vec) = self.map_of_moles_for_each_substance.get(*key) {
                    if row_idx < vec.len() {
                        row.push(Cell::new(&format!("{:.6}", vec[row_idx])));
                    } else {
                        row.push(Cell::new(""));
                    }
                } else {
                    row.push(Cell::new(""));
                }
            }
            table.add_row(Row::new(row));
        }

        table
    }

    /// Solve the equilibrium and return an error enum on failure
    pub fn solve(&mut self) -> Result<(), ReactionExtentError> {
        self.check_task()?;

        let (initial_guess, generated_initial_guess) = self.resolved_initial_guess()?;
        let generated_seed_to_publish = generated_initial_guess.then(|| initial_guess.clone());
        let candidate = self.solve_candidate_from_seed(initial_guess)?;
        if let Some(generated_seed) = generated_seed_to_publish {
            self.initial_guess = Some(generated_seed);
        }
        self.last_phase_control_report = None;
        self.publish_solve_candidate(candidate);
        Ok(())
    }

    /// Solves and validates one fixed formulation without publishing it.
    ///
    /// The supplied seed is explicit so an outer active-set algorithm cannot
    /// accidentally fall back to the seed from an earlier phase assemblage.
    pub(crate) fn solve_candidate_from_seed(
        &mut self,
        initial_guess: Vec<f64>,
    ) -> Result<EquilibriumSolveCandidate, ReactionExtentError> {
        self.check_task()?;
        if initial_guess.len() != self.n0.len()
            || initial_guess.iter().any(|value| !value.is_finite())
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_log_moles",
                message: format!(
                    "candidate seed must contain {} finite values, received {}",
                    self.n0.len(),
                    initial_guess.len()
                ),
            });
        }
        let policy = self
            .solver_settings
            .solver_policy
            .clone()
            .unwrap_or_else(|| {
                if self.has_rst_symbolic_context() {
                    SolverPolicy::rusted_scithe_default()
                } else {
                    SolverPolicy::legacy_default(self.solver_settings.solver)
                }
            });
        if policy.ordered_backends().is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "a solver policy must contain at least one backend".to_string(),
            });
        }

        let contract = self.build_solve_contract()?;
        info!("Initial guess: {:?}", &initial_guess);
        info!("n0 (species initial moles): {:?}", &self.n0);
        info!("Stoichiometric matrix:\n{}", &self.stoich_matrix);

        let try_step = contract.raw_residual(&initial_guess)?;
        if try_step.iter().any(|v| v.is_nan() || v.is_infinite()) {
            error!("Initial guess produces invalid residuals: {:?}", &try_step);
            return Err(ReactionExtentError::InvalidInitialResiduals(try_step));
        }

        let feasible = {
            move |xi: &[f64]| {
                compute_species_moles(xi)
                    .map(|moles| moles.iter().all(|&moles_i| moles_i >= -1e-12))
                    .unwrap_or(false)
            }
        };

        let validate_candidate = |candidate: &[f64]| {
            let raw_residual = contract.raw_residual(candidate)?;
            let acceptance_residual = contract.residual(candidate)?;
            validate_equilibrium_candidate(
                EquilibriumCandidateResiduals {
                    log_moles: candidate,
                    raw_residual: &raw_residual,
                    acceptance_residual: &acceptance_residual,
                },
                EquilibriumAcceptanceCriteria::new(
                    self.solver_settings.solver_params.tol,
                    10.0 * self.solver_settings.solver_params.tol.max(1e-8),
                    self.solver_settings.solver_params.tol,
                )?,
                &self.elem_composition,
                &self.elements_vector,
            )
        };

        let (sol, validation_report, solve_report) = self.solver_impl(
            initial_guess,
            contract.f.as_ref(),
            contract.jacobian(),
            &feasible,
            &validate_candidate,
            policy,
            contract.rst_problem.as_ref(),
        )?;

        let (moles, mole_table) = self.reconstructed_mole_state(&sol)?;
        let conditions = EquilibriumConditions::new(self.T, self.P, self.p0)?;
        let candidate_solution = EquilibriumSolution::new(
            sol.clone(),
            moles.clone(),
            conditions,
            validation_report.clone(),
        )?;
        let keq_validation_status = self.run_keq_cross_validation(&candidate_solution)?;
        Ok(EquilibriumSolveCandidate {
            log_moles: sol,
            moles,
            mole_table,
            validation_report,
            solve_report,
            keq_validation_status,
        })
    }

    /// Atomically publishes one candidate that has already passed every
    /// backend-independent validation gate.
    pub(crate) fn publish_solve_candidate(&mut self, candidate: EquilibriumSolveCandidate) {
        self.last_keq_validation_status = candidate.keq_validation_status;
        self.moles = candidate.moles;
        self.map_of_moles_for_each_substance = candidate.mole_table;
        self.solution = candidate.log_moles;
        self.last_validation_report = Some(candidate.validation_report);
        self.last_solve_report = Some(candidate.solve_report);
    }
    /// Runs the explicit backend policy and records every numerical outcome.
    pub(crate) fn solver_impl(
        &self,
        initial_guess: Vec<f64>,
        f: &dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
        j: Option<&dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
        feasible: &dyn Fn(&[f64]) -> bool,
        validate_candidate: &dyn Fn(
            &[f64],
        )
            -> Result<EquilibriumCandidateReport, ReactionExtentError>,
        policy: SolverPolicy,
        rst_problem: Option<&RstPreparedProblem>,
    ) -> Result<(Vec<f64>, EquilibriumCandidateReport, EquilibriumSolveReport), ReactionExtentError>
    {
        let backends = policy.ordered_backends();
        if backends.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "a solver policy must contain at least one backend".to_string(),
            });
        }

        let budget = self.solver_settings.solver_budget.unwrap_or_else(|| {
            SolverCascadeBudget::new(
                backends.len(),
                self.solver_settings.solver_params.max_iter,
                self.solver_settings
                    .solver_params
                    .max_iter
                    .saturating_mul(backends.len()),
            )
        });
        let backend_refs: Vec<&dyn EquilibriumNonlinearBackend> = backends
            .iter()
            .map(|backend| backend as &dyn EquilibriumNonlinearBackend)
            .collect();

        Self::solve_backend_cascade(
            &backend_refs,
            initial_guess,
            f,
            j,
            feasible,
            validate_candidate,
            policy,
            budget,
            &self.solver_settings.solver_params,
            &self.n0,
            &self.stoich_matrix,
            rst_problem,
        )
    }

    /// Runs the explicit backend policy and records every numerical outcome.
    pub(crate) fn solve_backend_cascade(
        backends: &[&dyn EquilibriumNonlinearBackend],
        initial_guess: Vec<f64>,
        f: &dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
        j: Option<&dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
        feasible: &dyn Fn(&[f64]) -> bool,
        validate_candidate: &dyn Fn(
            &[f64],
        )
            -> Result<EquilibriumCandidateReport, ReactionExtentError>,
        policy: SolverPolicy,
        budget: SolverCascadeBudget,
        solver_params: &SolverParams,
        initial_moles: &[f64],
        reactions: &DMatrix<f64>,
        rst_problem: Option<&RstPreparedProblem>,
    ) -> Result<(Vec<f64>, EquilibriumCandidateReport, EquilibriumSolveReport), ReactionExtentError>
    {
        if backends.is_empty() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_policy",
                message: "a solver policy must contain at least one backend".to_string(),
            });
        }

        if budget.max_attempts == 0
            || budget.max_iterations_per_attempt == 0
            || budget.max_total_iterations == 0
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "solver_budget",
                message: "attempt, per-attempt, and total iteration limits must be positive"
                    .to_string(),
            });
        }

        let mut attempts = Vec::with_capacity(backends.len());
        let mut started_attempts = 0usize;
        let mut remaining_iterations = budget.max_total_iterations;
        for &backend in backends {
            let backend_id = backend.backend();
            if started_attempts >= budget.max_attempts || remaining_iterations == 0 {
                let reason = if started_attempts >= budget.max_attempts {
                    format!(
                        "cascade attempt budget exhausted after {} started backend(s)",
                        budget.max_attempts
                    )
                } else {
                    "cascade total iteration budget exhausted".to_string()
                };
                attempts.push(SolverAttemptReport {
                    backend: backend_id,
                    outcome: SolverAttemptOutcome::Skipped { reason },
                    metrics: None,
                });
                continue;
            }
            let max_iterations = budget.max_iterations_per_attempt.min(remaining_iterations);
            started_attempts += 1;
            remaining_iterations -= max_iterations;

            // Every backend receives this same validated seed. A previous
            // rejected candidate is evidence for diagnostics, never an
            // implicit warm start.
            let result = backend.solve(BackendSolveRequest {
                initial_guess: initial_guess.clone(),
                residual: f,
                jacobian: j,
                feasible,
                params: solver_params,
                initial_moles,
                reactions,
                max_iterations,
                rst_problem,
            });

            match result {
                Ok(candidate) => match validate_candidate(&candidate.solution) {
                    Ok(validation) => {
                        attempts.push(SolverAttemptReport {
                            backend: backend_id,
                            outcome: SolverAttemptOutcome::Accepted,
                            metrics: candidate.metrics,
                        });
                        let report = EquilibriumSolveReport {
                            policy,
                            attempts,
                            accepted_backend: backend_id,
                        };
                        if report.accepted_after_fallback() {
                            info!("Solver accepted after fallback: {}", report);
                        }
                        return Ok((candidate.solution, validation, report));
                    }
                    Err(error) => attempts.push(SolverAttemptReport {
                        backend: backend_id,
                        outcome: SolverAttemptOutcome::RejectedCandidate {
                            reason: error.to_string(),
                        },
                        metrics: candidate.metrics,
                    }),
                },
                Err(e) => {
                    let Some(kind) = recoverable_backend_failure_kind(&e) else {
                        return if attempts.is_empty() {
                            Err(e)
                        } else {
                            Err(ReactionExtentError::CascadeAborted {
                                attempts,
                                cause: Box::new(e),
                            })
                        };
                    };
                    attempts.push(SolverAttemptReport {
                        backend: backend_id,
                        outcome: SolverAttemptOutcome::Failed {
                            kind,
                            // Attempt reports are part of the public solve
                            // contract, so keep their wording independent of
                            // a private Debug implementation.
                            reason: e.to_string(),
                        },
                        metrics: None,
                    });
                }
            }
        }

        Err(ReactionExtentError::AllBackendsFailed { attempts })
    }

    /// Converts log-mole solution to actual mole numbers
    ///
    /// Transforms the solver output (ln(n_i)) back to mole numbers (n_i = exp(ln(n_i))).
    pub fn compute_species_moles(&mut self, sol: Vec<f64>) -> Result<(), ReactionExtentError> {
        let (moles, map_of_moles_for_each_substance) = self.reconstructed_mole_state(&sol)?;
        self.moles = moles;
        self.map_of_moles_for_each_substance = map_of_moles_for_each_substance;
        Ok(())
    }

    /// Reconstructs the complete public mole view without mutating solver
    /// state. Both regular solve publication and the compatibility mutator use
    /// this one validation boundary.
    pub(crate) fn reconstructed_mole_state(
        &self,
        sol: &[f64],
    ) -> Result<(Vec<f64>, HashMap<String, Vec<f64>>), ReactionExtentError> {
        let moles = compute_species_moles(&sol)?;
        if moles.len() != self.subs_data.substances.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "solution has {} mole values for {} substances",
                moles.len(),
                self.subs_data.substances.len(),
            )));
        }
        let map_of_moles_for_each_substance: HashMap<String, Vec<f64>> = self
            .subs_data
            .substances
            .clone()
            .iter()
            .zip(moles.clone())
            .map(|(subs, n)| (subs.clone(), vec![n]))
            .collect();
        if map_of_moles_for_each_substance.len() != moles.len() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "species",
                message: "species names must be unique when publishing mole results".to_string(),
            });
        }
        Ok((moles, map_of_moles_for_each_substance))
    }

    /// Validates problem setup for dimensional consistency
    ///
    /// Ensures that initial moles, substances, and initial guess vectors
    /// all have consistent dimensions.
    pub(crate) fn check_task(&self) -> Result<(), ReactionExtentError> {
        self.solver_settings.validate()?;
        let n_subs = self.subs_data.substances.len();
        let mut seen_species = HashSet::with_capacity(n_subs);
        if self
            .subs_data
            .substances
            .iter()
            .any(|species| !seen_species.insert(species))
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "species",
                message: "species names must be unique".to_string(),
            });
        }
        if self.n0.len() != n_subs {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "initial moles number of substances {} doesn't
            match number of substances {}",
                self.n0.len(),
                n_subs
            )));
        }
        if self.gibbs.len() != n_subs {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "Gibbs vector has {} entries but there are {n_subs} substances",
                self.gibbs.len(),
            )));
        }
        if self
            .n0
            .iter()
            .any(|moles| !moles.is_finite() || *moles < 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "initial_moles",
                message: "initial moles must be finite and non-negative".to_string(),
            });
        }
        if let Some(initial_guess) = self.initial_guess.as_ref() {
            if initial_guess.len() != n_subs {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "initial guess for solution {} doesn't match number of substances {n_subs}",
                    initial_guess.len(),
                )));
            }
            if initial_guess.iter().any(|value| !value.is_finite()) {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "initial_log_moles",
                    message: "initial log-mole values must be finite".to_string(),
                });
            }
        }
        Ok(())
    }
}

/// A fallback policy can recover from numerical failures, but it must not
/// hide malformed input or invalid solver configuration behind later
/// backends. Candidate validation is handled separately and is retryable by
/// design because a later method can produce a better iterate.
#[cfg(test)]
fn is_recoverable_backend_failure(error: &ReactionExtentError) -> bool {
    recoverable_backend_failure_kind(error).is_some()
}

/// Classifies only failures for which another numerical backend may still
/// produce an accepted candidate. Invalid domain input deliberately maps to
/// `None` so a fallback can never hide a malformed equilibrium problem.
pub(crate) fn recoverable_backend_failure_kind(
    error: &ReactionExtentError,
) -> Option<SolverAttemptFailureKind> {
    if !error.is_retryable_backend_failure() {
        return None;
    }
    match error {
        ReactionExtentError::SolveError(_) => Some(SolverAttemptFailureKind::Solver),
        ReactionExtentError::BackendFailure { .. } => Some(SolverAttemptFailureKind::Backend),
        ReactionExtentError::ResidualEvaluation(_) => {
            Some(SolverAttemptFailureKind::ResidualEvaluation)
        }
        ReactionExtentError::JacobianEvaluation(_) => {
            Some(SolverAttemptFailureKind::JacobianEvaluation)
        }
        _ => None,
    }
}

#[cfg(test)]
mod retry_policy_tests {
    use super::{is_recoverable_backend_failure, recoverable_backend_failure_kind};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        BackendFailureKind, ReactionExtentError, SolveError,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptFailureKind;

    #[test]
    fn retries_only_numerical_backend_failures() {
        assert!(is_recoverable_backend_failure(
            &ReactionExtentError::SolveError(SolveError::MaxIterations,)
        ));
        assert!(is_recoverable_backend_failure(
            &ReactionExtentError::JacobianEvaluation("singular numerical step".to_string()),
        ));
        assert!(is_recoverable_backend_failure(
            &ReactionExtentError::BackendFailure {
                backend: "rst_levenberg_marquardt".to_string(),
                kind: BackendFailureKind::SingularJacobian,
                message: "ill-conditioned system".to_string(),
            },
        ));
        assert!(!is_recoverable_backend_failure(
            &ReactionExtentError::InvalidProblem {
                field: "rst_solve_options",
                message: "max_iterations must be greater than zero".to_string(),
            },
        ));
        assert!(!is_recoverable_backend_failure(
            &ReactionExtentError::DimensionMismatch("species count".to_string()),
        ));
        assert!(
            ReactionExtentError::InvalidProblem {
                field: "solver_options",
                message: "invalid settings".to_string(),
            }
            .is_non_retryable_input_error()
        );
        assert_eq!(
            recoverable_backend_failure_kind(&ReactionExtentError::ResidualEvaluation(
                "temporary residual failure".to_string(),
            )),
            Some(SolverAttemptFailureKind::ResidualEvaluation),
        );
    }
}

#[cfg(test)]
mod solver_budget_tests {
    use super::{EquilibriumLogMoles, temperature_failure};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        SolverAttemptFailureKind, SolverAttemptOutcome, SolverAttemptReport, SolverBackend,
        SolverCascadeBudget, SolverPolicy,
    };
    use nalgebra::DMatrix;
    use std::cell::RefCell;

    #[test]
    fn global_iteration_budget_skips_later_backends_in_declared_order() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.set_solver_budget(SolverCascadeBudget::new(2, 3, 1));
        let forced_failure = |_values: &[f64]| {
            Err(ReactionExtentError::ResidualEvaluation(
                "forced numerical failure".to_string(),
            ))
        };
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        let accept = |_values: &[f64]| {
            Err(ReactionExtentError::InvalidCandidate {
                field: "test_acceptance_gate",
                message: "candidate should not be reached".to_string(),
            })
        };

        let error = solver
            .solver_impl(
                vec![0.0],
                &forced_failure,
                Some(&jacobian),
                &feasible,
                &accept,
                SolverPolicy::Cascade(vec![
                    SolverBackend::Legacy(Solvers::LM),
                    SolverBackend::Legacy(Solvers::NR),
                ]),
                None,
            )
            .unwrap_err();
        let attempts = match error {
            ReactionExtentError::AllBackendsFailed { attempts } => attempts,
            other => panic!("expected a cascade report, got {other:?}"),
        };

        assert_eq!(attempts.len(), 2);
        let SolverAttemptOutcome::Failed { kind, reason } = &attempts[0].outcome else {
            panic!(
                "expected the first backend to fail, got {:?}",
                attempts[0].outcome
            );
        };
        assert_eq!(*kind, SolverAttemptFailureKind::Solver);
        assert_eq!(
            reason,
            "equilibrium solver failed: nonlinear solver evaluation failed: residual evaluation failed: equilibrium residual evaluation failed: forced numerical failure"
        );
        assert!(matches!(
            attempts[1].outcome,
            SolverAttemptOutcome::Skipped { .. }
        ));
    }

    #[test]
    fn every_started_backend_receives_the_original_validated_seed() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.set_solver_budget(SolverCascadeBudget::new(2, 1, 2));
        let observed_seeds = RefCell::new(Vec::new());
        let forced_failure = |values: &[f64]| {
            observed_seeds.borrow_mut().push(values.to_vec());
            Err(ReactionExtentError::ResidualEvaluation(
                "forced numerical failure".to_string(),
            ))
        };
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        let accept = |_values: &[f64]| {
            Err(ReactionExtentError::InvalidCandidate {
                field: "test_acceptance_gate",
                message: "candidate should not be reached".to_string(),
            })
        };

        let _ = solver.solver_impl(
            vec![-3.0],
            &forced_failure,
            Some(&jacobian),
            &feasible,
            &accept,
            SolverPolicy::Cascade(vec![
                SolverBackend::Legacy(Solvers::LM),
                SolverBackend::Legacy(Solvers::NR),
            ]),
            None,
        );

        assert_eq!(observed_seeds.into_inner(), vec![vec![-3.0], vec![-3.0]]);
    }

    #[test]
    fn non_retryable_cascade_abort_retains_earlier_attempts_and_root_cause() {
        let solver = EquilibriumLogMoles::empty();
        let recoverable_failure = |_values: &[f64]| {
            Err(ReactionExtentError::ResidualEvaluation(
                "forced numerical failure".to_string(),
            ))
        };
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        let accept = |_values: &[f64]| {
            Err(ReactionExtentError::InvalidCandidate {
                field: "test_acceptance_gate",
                message: "candidate should not be reached".to_string(),
            })
        };

        let error = solver
            .solver_impl(
                vec![0.0],
                &recoverable_failure,
                Some(&jacobian),
                &feasible,
                &accept,
                SolverPolicy::Cascade(vec![
                    SolverBackend::Legacy(Solvers::LM),
                    SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
                ]),
                None,
            )
            .unwrap_err();

        let ReactionExtentError::CascadeAborted { attempts, cause } = error else {
            panic!("expected the cascade to retain its abort trace");
        };
        assert_eq!(attempts.len(), 1);
        assert!(matches!(
            attempts[0].outcome,
            SolverAttemptOutcome::Failed { .. }
        ));
        assert!(matches!(
            *cause,
            ReactionExtentError::InvalidProblem {
                field: "rst_symbolic_problem",
                ..
            }
        ));
    }

    #[test]
    fn temperature_failure_retains_the_cascade_trace() {
        let attempts = vec![SolverAttemptReport {
            backend: SolverBackend::Legacy(Solvers::LM),
            outcome: SolverAttemptOutcome::Skipped {
                reason: "iteration budget exhausted".to_string(),
            },
            metrics: None,
        }];
        let failure = temperature_failure(
            1200.0,
            &ReactionExtentError::AllBackendsFailed {
                attempts: attempts.clone(),
            },
        );

        assert_eq!(failure.temperature, 1200.0);
        assert_eq!(failure.attempts, attempts);
        assert_eq!(
            failure.message,
            "all equilibrium solver backends failed after 1 attempt(s)"
        );

        let aborted = temperature_failure(
            1300.0,
            &ReactionExtentError::CascadeAborted {
                attempts: attempts.clone(),
                cause: Box::new(ReactionExtentError::InvalidProblem {
                    field: "rst_symbolic_problem",
                    message: "symbolic preparation failed".to_string(),
                }),
            },
        );
        assert_eq!(aborted.attempts, attempts);
        assert!(aborted.message.contains("cascade aborted"));
    }
}

#[cfg(test)]
mod solver_cascade_story_tests {
    use super::{EquilibriumLogMoles, SolverParams};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_backend_adapter::{
        BackendSolveRequest, BackendSolveResult, EquilibriumNonlinearBackend,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RustedSciTheSolver;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
        EquilibriumSolveReport, SolverAttemptFailureKind, SolverAttemptOutcome, SolverBackend,
        SolverCascadeBudget, SolverPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_validation::EquilibriumCandidateReport;
    use nalgebra::DMatrix;

    #[derive(Clone)]
    enum FakeBackendBehavior {
        RecoverableFailure,
        Success { solution: Vec<f64> },
    }

    struct FakeBackend {
        backend: SolverBackend,
        behavior: FakeBackendBehavior,
    }

    impl EquilibriumNonlinearBackend for FakeBackend {
        fn backend(&self) -> SolverBackend {
            self.backend
        }

        fn solve(
            &self,
            request: BackendSolveRequest<'_>,
        ) -> Result<BackendSolveResult, ReactionExtentError> {
            match &self.behavior {
                FakeBackendBehavior::RecoverableFailure => {
                    Err(ReactionExtentError::ResidualEvaluation(
                        "forced recoverable failure".to_string(),
                    ))
                }
                FakeBackendBehavior::Success { solution } => {
                    let _ = request;
                    Ok(BackendSolveResult {
                        solution: solution.clone(),
                        metrics: Some(
                            crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverAttemptMetrics {
                                termination: crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverTermination::Converged,
                                backend_converged: true,
                                iterations: 1,
                                residual_evaluations: 1,
                                jacobian_evaluations: 1,
                                linear_solves: 1,
                                elapsed_millis: 0,
                            },
                        ),
                    })
                }
            }
        }
    }

    fn accepted_candidate_report() -> EquilibriumCandidateReport {
        EquilibriumCandidateReport {
            residual_l2_norm: 0.0,
            residual_rms: 0.0,
            max_abs_residual: 0.0,
            raw_residual_l2_norm: 0.0,
            raw_residual_rms: 0.0,
            raw_max_abs_residual: 0.0,
            max_abs_element_balance_error: 0.0,
            reaction_affinity_l2_norm: 0.0,
            max_abs_reaction_affinity: 0.0,
            min_moles: 1.0,
        }
    }

    fn solve_with_fake_backends(
        backends: Vec<&dyn EquilibriumNonlinearBackend>,
        validate_candidate: impl Fn(&[f64]) -> Result<EquilibriumCandidateReport, ReactionExtentError>,
        budget: SolverCascadeBudget,
    ) -> Result<(Vec<f64>, EquilibriumCandidateReport, EquilibriumSolveReport), ReactionExtentError>
    {
        let residual = |values: &[f64]| Ok(values.to_vec());
        let jacobian = |_values: &[f64]| Ok(DMatrix::identity(1, 1));
        let feasible = |_values: &[f64]| true;
        let policy =
            SolverPolicy::Cascade(backends.iter().map(|backend| backend.backend()).collect());

        EquilibriumLogMoles::solve_backend_cascade(
            &backends,
            vec![0.0],
            &residual,
            Some(&jacobian),
            &feasible,
            &validate_candidate,
            policy,
            budget,
            &SolverParams::default(),
            &[1.0],
            &DMatrix::from_row_slice(1, 1, &[1.0]),
            None,
        )
    }

    #[test]
    fn first_backend_fails_recoverably_second_succeeds_and_is_accepted() {
        let fail = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::LM),
            behavior: FakeBackendBehavior::RecoverableFailure,
        };
        let success = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::NR),
            behavior: FakeBackendBehavior::Success {
                solution: vec![2.0],
            },
        };

        let result = solve_with_fake_backends(
            vec![&fail, &success],
            |_candidate| Ok(accepted_candidate_report()),
            SolverCascadeBudget::new(2, 3, 6),
        )
        .expect("the second backend should be accepted");

        assert_eq!(result.0, vec![2.0]);
        assert_eq!(result.2.attempts.len(), 2);
        assert!(matches!(
            result.2.attempts[0].outcome,
            SolverAttemptOutcome::Failed {
                kind: SolverAttemptFailureKind::ResidualEvaluation,
                ..
            }
        ));
        assert_eq!(result.2.attempts[1].outcome, SolverAttemptOutcome::Accepted);
        assert_eq!(
            result.2.accepted_backend,
            SolverBackend::Legacy(super::Solvers::NR)
        );
    }

    #[test]
    fn first_backend_returns_a_rejected_candidate_and_the_next_one_succeeds() {
        let first = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::LM),
            behavior: FakeBackendBehavior::Success {
                solution: vec![1.0],
            },
        };
        let second = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::NR),
            behavior: FakeBackendBehavior::Success {
                solution: vec![2.0],
            },
        };

        let result = solve_with_fake_backends(
            vec![&first, &second],
            |candidate| {
                if candidate[0] == 1.0 {
                    Err(ReactionExtentError::InvalidCandidate {
                        field: "candidate_acceptance",
                        message: "first backend candidate is intentionally rejected".to_string(),
                    })
                } else {
                    Ok(accepted_candidate_report())
                }
            },
            SolverCascadeBudget::new(2, 3, 6),
        )
        .expect("the second backend should be accepted after rejection");

        assert_eq!(result.0, vec![2.0]);
        assert_eq!(result.2.attempts.len(), 2);
        assert!(matches!(
            result.2.attempts[0].outcome,
            SolverAttemptOutcome::RejectedCandidate { .. }
        ));
        assert_eq!(result.2.attempts[1].outcome, SolverAttemptOutcome::Accepted);
    }

    #[test]
    fn every_backend_fails_and_the_report_keeps_the_declared_order() {
        let first = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::LM),
            behavior: FakeBackendBehavior::RecoverableFailure,
        };
        let second = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::NR),
            behavior: FakeBackendBehavior::RecoverableFailure,
        };

        let error = solve_with_fake_backends(
            vec![&first, &second],
            |_candidate| Ok(accepted_candidate_report()),
            SolverCascadeBudget::new(2, 3, 6),
        )
        .unwrap_err();

        let ReactionExtentError::AllBackendsFailed { attempts } = error else {
            panic!("expected an exhaustive cascade failure");
        };
        assert_eq!(attempts.len(), 2);
        assert_eq!(
            attempts[0].backend,
            SolverBackend::Legacy(super::Solvers::LM)
        );
        assert_eq!(
            attempts[1].backend,
            SolverBackend::Legacy(super::Solvers::NR)
        );
        assert!(matches!(
            attempts[0].outcome,
            SolverAttemptOutcome::Failed { .. }
        ));
        assert!(matches!(
            attempts[1].outcome,
            SolverAttemptOutcome::Failed { .. }
        ));
    }

    #[test]
    fn global_budget_exhaustion_stops_later_backends_deterministically() {
        let first = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::LM),
            behavior: FakeBackendBehavior::RecoverableFailure,
        };
        let second = FakeBackend {
            backend: SolverBackend::Legacy(super::Solvers::NR),
            behavior: FakeBackendBehavior::Success {
                solution: vec![2.0],
            },
        };
        let third = FakeBackend {
            backend: SolverBackend::RustedSciThe(RustedSciTheSolver::LevenbergMarquardt),
            behavior: FakeBackendBehavior::Success {
                solution: vec![3.0],
            },
        };

        let error = solve_with_fake_backends(
            vec![&first, &second, &third],
            |_candidate| Ok(accepted_candidate_report()),
            SolverCascadeBudget::new(1, 2, 2),
        )
        .unwrap_err();

        let ReactionExtentError::AllBackendsFailed { attempts } = error else {
            panic!("expected a budget-limited cascade failure");
        };
        assert_eq!(attempts.len(), 3);
        assert!(matches!(
            attempts[0].outcome,
            SolverAttemptOutcome::Failed { .. }
        ));
        assert!(matches!(
            attempts[1].outcome,
            SolverAttemptOutcome::Skipped { .. }
        ));
        assert!(matches!(
            attempts[2].outcome,
            SolverAttemptOutcome::Skipped { .. }
        ));
    }
}

#[cfg(test)]
mod continuation_seed_policy_tests {
    use super::{
        ContinuationSeedPolicy, EquilibriumLogMoles, Solvers, continuation_seed_for_point,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
        DEFAULT_TRACE_MOLE_FLOOR, TraceSpeciesSeedPolicy,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::gas_solver;

    #[test]
    fn default_policy_continues_from_the_last_accepted_point() {
        assert_eq!(
            EquilibriumLogMoles::empty()
                .solver_settings
                .continuation_seed_policy,
            ContinuationSeedPolicy::PreviousAccepted
        );
        assert_eq!(
            continuation_seed_for_point(
                ContinuationSeedPolicy::PreviousAccepted,
                &Some(vec![-1.0]),
                &Some(vec![-4.0]),
            ),
            Some(vec![-4.0])
        );
    }

    #[test]
    fn independent_policy_preserves_the_configured_seed() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.set_continuation_seed_policy(ContinuationSeedPolicy::IndependentPerPoint);

        assert_eq!(
            continuation_seed_for_point(
                solver.solver_settings.continuation_seed_policy,
                &Some(vec![-1.0]),
                &Some(vec![-4.0]),
            ),
            Some(vec![-1.0])
        );
    }

    #[test]
    fn solver_settings_use_the_absolute_trace_seed_by_default() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.n0 = vec![1.0, 0.0];

        let (seed, generated) = solver.resolved_initial_guess().unwrap();
        assert!(generated);
        assert_eq!(seed[0], 0.0);
        assert_eq!(seed[1], DEFAULT_TRACE_MOLE_FLOOR.ln());
    }

    #[test]
    fn solver_settings_can_generate_a_relative_trace_seed_explicitly() {
        let mut solver = EquilibriumLogMoles::empty();
        solver.n0 = vec![4.0, 0.0];
        solver.solver_settings.trace_seed_policy =
            TraceSpeciesSeedPolicy::RelativeToLargestInitialMole {
                fraction: 1e-6,
                minimum_floor: 1e-30,
            };
        solver.solver_settings.solver = Solvers::LM;

        let (seed, generated) = solver.resolved_initial_guess().unwrap();
        assert!(generated);
        assert_eq!(seed[0], 4.0_f64.ln());
        assert!((seed[1] - (4.0e-6_f64).ln()).abs() < 1e-12);
    }

    #[test]
    fn gibbs_closure_builders_preserve_species_order_for_both_rc_and_arc_wrappers() {
        let mut solver = gas_solver(
            vec!["O2".to_string(), "O".to_string()],
            1000.0,
            101325.0,
            Solvers::LM,
            None,
            false,
        )
        .unwrap();

        let serial = solver.build_gibbs_functions().unwrap();
        let parallel = solver.build_parallel_gibbs_functions().unwrap();

        assert_eq!(serial.len(), 2);
        assert_eq!(parallel.len(), 2);
        for temperature in [300.0, 1200.0, 2400.0] {
            for index in 0..serial.len() {
                let rc_value = serial[index](temperature);
                let arc_value = parallel[index](temperature);
                assert!(
                    (rc_value - arc_value).abs() <= 1e-12,
                    "builder mismatch for species {index} at T={temperature}"
                );
            }
        }
    }
}

#[cfg(test)]
mod activity_contract_regression_tests {
    #![allow(deprecated)]
    use super::{
        GibbsFn, Phase, PhaseKind, equilibrium_logmole_residual2,
        evaluate_equilibrium_logmole_residual, reaction_phase_stoichiometry,
    };
    use nalgebra::DMatrix;
    use std::rc::Rc;

    #[test]
    fn deprecated_residual_path_matches_the_canonical_activity_contract() {
        let reactions = DMatrix::from_row_slice(2, 1, &[-1.0, 1.0]);
        let elements = DMatrix::from_row_slice(2, 1, &[1.0, 1.0]);
        let element_totals = vec![2.0];
        let gibbs: Vec<GibbsFn> = vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0)];
        let phases = vec![
            Phase {
                kind: PhaseKind::IdealGas,
                species: vec![0],
            },
            Phase {
                kind: PhaseKind::IdealSolution,
                species: vec![1],
            },
        ];
        let species_phase = vec![0, 1];
        let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
        let y = vec![0.0, 0.0];

        let canonical = evaluate_equilibrium_logmole_residual(
            &y,
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            500.0,
            202_650.0,
            101_325.0,
            &species_phase,
            &delta_n,
        )
        .unwrap();
        let legacy = equilibrium_logmole_residual2(
            reactions,
            elements,
            element_totals,
            gibbs,
            phases,
            500.0,
            202_650.0,
            101_325.0,
            1e-30,
            1e-30,
        )
        .unwrap();
        let legacy = legacy(&y).unwrap();

        assert_eq!(canonical.len(), legacy.len());
        for (lhs, rhs) in canonical.iter().zip(legacy.iter()) {
            assert!((lhs - rhs).abs() <= 1e-12);
        }
    }
}

/// Reconstructs finite, strictly positive physical mole values from log-moles.
///
/// The logarithmic formulation intentionally represents trace species, but it
/// must not silently turn an extreme iterate into an underflowed zero or an
/// overflowed infinity. Callers decide whether this is a candidate, residual,
/// or Jacobian failure; the reconstruction contract itself stays uniform.
pub fn compute_species_moles(log_moles: &[f64]) -> Result<Vec<f64>, ReactionExtentError> {
    let mut moles = Vec::with_capacity(log_moles.len());
    for (index, &log_moles_i) in log_moles.iter().enumerate() {
        if !log_moles_i.is_finite() {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "log_moles",
                message: format!("log-mole coordinate {index} is non-finite"),
            });
        }
        let moles_i = log_moles_i.exp();
        if !moles_i.is_finite() || moles_i <= 0.0 {
            return Err(ReactionExtentError::InvalidCandidate {
                field: "mole_amount",
                message: format!(
                    "log-mole coordinate {index} = {log_moles_i:e} reconstructs to an invalid mole amount"
                ),
            });
        }
        moles.push(moles_i);
    }
    Ok(moles)
}

/// Builds the column-oriented sweep table without changing any solver field.
///
/// Sweep entry points use this as a final validation step before publishing a
/// bundle of rows, failures, snapshots, and the derived table together.
fn build_mole_series_table(
    substances: &[String],
    rows: &[(f64, Vec<f64>)],
) -> Result<HashMap<String, Vec<f64>>, ReactionExtentError> {
    let mut unique_species = HashSet::with_capacity(substances.len());
    if substances
        .iter()
        .any(|substance| !unique_species.insert(substance))
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "species",
            message: "species names must be unique when publishing sweep results".to_string(),
        });
    }
    for (temperature, moles) in rows {
        if moles.len() != substances.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "temperature {temperature} produced {} mole values for {} substances",
                moles.len(),
                substances.len(),
            )));
        }
    }

    let mut published = HashMap::with_capacity(substances.len() + 1);
    published.insert(
        "T".to_string(),
        rows.iter().map(|(temperature, _)| *temperature).collect(),
    );
    for (index, substance) in substances.iter().enumerate() {
        published.insert(
            substance.clone(),
            rows.iter().map(|(_, moles)| moles[index]).collect(),
        );
    }
    Ok(published)
}

/// Computes total element amounts from species composition.
///
/// The check is deliberately performed before `nalgebra` multiplication: a
/// mismatched public matrix/vector pair must become a typed domain error, not
/// a linear-algebra assertion failure.
pub fn compute_element_totals(
    a: &DMatrix<f64>,
    n0: &[f64],
) -> Result<DVector<f64>, ReactionExtentError> {
    if a.nrows() != n0.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "element composition has {} species rows but initial moles have {} entries",
            a.nrows(),
            n0.len(),
        )));
    }
    if a.iter().any(|value| !value.is_finite()) || n0.iter().any(|value| !value.is_finite()) {
        return Err(ReactionExtentError::InvalidProblem {
            field: "element_composition_or_initial_moles",
            message: "element composition and initial moles must be finite".to_string(),
        });
    }
    Ok(&a.transpose() * DVector::from_column_slice(n0))
}

/// Calculates standard Gibbs free energy change for each reaction
///
/// Computes ΔG°_rxn = Σ ν_i * G°_i for each independent reaction.
pub fn reaction_standard_gibbs(stoich: &DMatrix<f64>, gibbs: &[GibbsFn], T: f64) -> Vec<f64> {
    let m = stoich.nrows();
    let r = stoich.ncols();

    let mut dg0 = vec![0.0; r];

    for k in 0..r {
        let mut sum = 0.0;
        for i in 0..m {
            sum += stoich[(i, k)] * gibbs[i](T);
        }
        dg0[k] = sum;
    }
    dg0
}

/// Computes row scaling factors for the canonical log-mole residual.
///
/// The returned vector has one entry per residual row:
/// - the first `r` entries correspond to reaction-equilibrium equations and
///   are measured in the same dimensionless log-space as the residual itself;
/// - the remaining `e` entries correspond to element-balance equations and are
///   measured in mole units.
///
/// Scaling is a numerical conditioning aid only. It does not change the
/// physical variables, and it must be applied to residual rows and Jacobian
/// rows consistently.
pub fn equilibrium_scaling(
    stoich: &DMatrix<f64>,   // m × r
    elements: &DMatrix<f64>, // m × E
    gibbs: &[GibbsFn],       // m
    element_totals: &[f64],  // E
    T: f64,
) -> Result<Vec<f64>, ReactionExtentError> {
    let m = stoich.nrows();
    let r = stoich.ncols();
    let e = elements.ncols();
    if !T.is_finite() || T <= 0.0 {
        return Err(ReactionExtentError::InvalidConditions {
            parameter: "temperature",
            value: T,
        });
    }
    if elements.nrows() != m {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "stoichiometric matrix has {m} species rows but element matrix has {}",
            elements.nrows(),
        )));
    }
    if gibbs.len() != m {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "scaling received {} Gibbs functions for {m} species",
            gibbs.len(),
        )));
    }
    if element_totals.len() != e {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "scaling received {} element totals for {e} elements",
            element_totals.len(),
        )));
    }

    let rt = 8.314462618 * T;

    // --- reaction scaling ---
    let dg0 = reaction_standard_gibbs(stoich, gibbs, T);
    let mut scale = vec![0.0; r + e];

    for k in 0..r {
        let mut nu_norm = 0.0;
        for i in 0..m {
            let nu = stoich[(i, k)];
            nu_norm += nu * nu;
        }
        nu_norm = nu_norm.sqrt();

        scale[k] = dg0[k].abs().max(rt * nu_norm).max(10.0); // dimensionless, log-scale
    }

    // --- element balance scaling ---
    for el in 0..e {
        scale[r + el] = element_totals[el].abs().max(10.0);
    }

    Ok(scale)
}

/// Legacy spelling retained only for fixture literals while callers migrate.
///
/// This is a type alias, not a second model. All solver production paths store
/// and consume [`PhaseActivityModel`] directly.
pub type PhaseKind = PhaseActivityModel;

/// Represents a thermodynamic phase containing multiple species
#[derive(Debug, Clone)]
pub struct Phase {
    /// Canonical activity law (ideal gas or ideal solution).
    pub kind: PhaseActivityModel,
    /// Indices of species belonging to this phase
    pub species: Vec<usize>,
}

/// Function type for Gibbs free energy evaluation
///
/// Takes temperature in K and returns standard Gibbs free energy in J/mol.
pub type GibbsFn = Rc<dyn Fn(f64) -> f64>;

/// Maps species indices to their corresponding phase indices
///
/// Validates that each species belongs to exactly one phase and returns
/// a vector where element i gives the phase index for species i.
pub fn species_to_phase_map(
    phases: &[Phase],
    num_species: usize,
) -> Result<Vec<usize>, ReactionExtentError> {
    let mut map = vec![None; num_species];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            if i >= num_species {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "phase {j} references species index {i}, but the problem has {num_species} species"
                )));
            }
            if map[i].is_some() {
                error!(
                    "Species {} appears in more than one phase. Duplicate species per phase instead (e.g. O2_g, O2_l).",
                    i
                );
                return Err(ReactionExtentError::DuplicateSpecies(i));
            }
            // Assign species i to phase j
            map[i] = Some(j);
        }
    }

    let mut out = Vec::with_capacity(num_species);
    for x in map.into_iter() {
        match x {
            Some(v) => out.push(v),
            None => {
                error!("Species not assigned to any phase");
                return Err(ReactionExtentError::SpeciesNotAssigned { index: out.len() });
            }
        }
    }
    Ok(out)
}

/// Computes change in mole numbers per phase for each reaction
///
/// For each reaction k and phase j, calculates Δn[k][j] = sum of stoichiometric
/// coefficients for all species in phase j participating in reaction k.
pub fn reaction_phase_stoichiometry(
    reactions: &DMatrix<f64>, // m × r
    phases: &[Phase],
) -> Vec<Vec<f64>> {
    let r = reactions.ncols();

    let mut delta_n = vec![vec![0.0; phases.len()]; r];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            for k in 0..r {
                delta_n[k][j] += reactions[(i, k)];
            }
        }
    }

    delta_n
}

// The threshold-based residual/Jacobian pair below is retained temporarily as
// a deprecated compatibility surface. The canonical solver never selects it;
// new code must use the validated log-mole formulation without activity masks.
#[inline]
fn species_active(n: f64, eps: f64) -> bool {
    n > eps
}

#[inline]
fn phase_active(n: f64, eps: f64) -> bool {
    n > eps
}

/// Validates the square log-mole equilibrium system before numerical closures
/// capture its data. The formulation has one equation per species: independent
/// reaction equations plus elemental-balance equations must equal `m`.
pub(crate) fn validate_logmole_system_dimensions(
    reactions: &DMatrix<f64>,
    elements: &DMatrix<f64>,
    species_phase: &[usize],
    phase_count: usize,
    delta_n: &[Vec<f64>],
) -> Result<(), ReactionExtentError> {
    let species_count = reactions.nrows();
    let reaction_count = reactions.ncols();
    let element_count = elements.ncols();

    if elements.nrows() != species_count {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "element matrix has {} species rows, expected {species_count}",
            elements.nrows()
        )));
    }
    if reaction_count + element_count != species_count {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "log-mole formulation has {reaction_count} reaction and {element_count} element equations for {species_count} species"
        )));
    }
    if phase_count == 0 || species_phase.len() != species_count {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "species-to-phase map has {} entries for {species_count} species and {phase_count} phases",
            species_phase.len()
        )));
    }
    if species_phase.iter().any(|&phase| phase >= phase_count) {
        return Err(ReactionExtentError::DimensionMismatch(
            "species-to-phase map references a missing phase".to_string(),
        ));
    }
    if delta_n.len() != reaction_count
        || delta_n
            .iter()
            .any(|phase_row| phase_row.len() != phase_count)
    {
        return Err(ReactionExtentError::DimensionMismatch(
            "reaction-to-phase stoichiometry dimensions do not match the problem".to_string(),
        ));
    }

    Ok(())
}

/// Validates the three scalar conditions consumed directly by the residual.
///
/// This is deliberately separate from matrix validation so every caller gets
/// the exact offending parameter and value instead of an ambiguous combined
/// "temperature or pressure" diagnostic.
pub(crate) fn validate_residual_conditions(
    temperature: f64,
    pressure: f64,
    reference_pressure: f64,
) -> Result<(), ReactionExtentError> {
    for (parameter, value) in [
        ("temperature", temperature),
        ("pressure", pressure),
        ("reference_pressure", reference_pressure),
    ] {
        if !value.is_finite() || value <= 0.0 {
            return Err(ReactionExtentError::InvalidConditions { parameter, value });
        }
    }
    Ok(())
}

/// Generates a legacy residual closure for the canonical log-mole system.
///
/// The residual has two blocks:
/// - reaction equations, written in dimensionless log-space;
/// - element-balance equations, written in mole units.
///
/// This closure factory exists for compatibility with the temporary legacy
/// backend adapter. New code should prefer the pure
/// [`evaluate_equilibrium_logmole_residual`] entry point on a prepared problem.
///
pub fn equilibrium_logmole_residual(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs: Vec<GibbsFn>,      // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    species_phase: Vec<usize>,
    _subs_eps: f64,
    __phase_eps: f64,
) -> Result<Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>, ReactionExtentError> {
    let m = reactions.nrows();
    let e = elements.ncols();

    if gibbs.len() != m || element_totals.len() != e {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "log-mole residual has {m} species, {} Gibbs functions, and {} element totals for {e} elements",
            gibbs.len(),
            element_totals.len(),
        )));
    }

    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    validate_logmole_system_dimensions(
        &reactions,
        &elements,
        &species_phase,
        phases.len(),
        &delta_n,
    )?;
    validate_residual_conditions(temperature, pressure, p0)?;
    Ok(Box::new(move |log_moles| {
        evaluate_equilibrium_logmole_residual(
            log_moles,
            &reactions,
            &elements,
            &element_totals,
            &gibbs,
            &phases,
            temperature,
            pressure,
            p0,
            &species_phase,
            &delta_n,
        )
    }))
}

/// Evaluates the canonical log-mole residual for one prepared state.
///
/// Residual layout:
/// - rows `0..r` are reaction-equilibrium equations in dimensionless log
///   space;
/// - rows `r..r+e` are element-balance equations in mole units.
///
/// The row ordering is deterministic and must match the prepared problem's
/// reaction basis and element layout exactly. This function is the canonical
/// implementation; legacy closure adapters delegate here rather than keeping a
/// second formula alive.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_equilibrium_logmole_residual(
    log_moles: &[f64],
    reactions: &DMatrix<f64>,
    elements: &DMatrix<f64>,
    element_totals: &[f64],
    gibbs: &[GibbsFn],
    phases: &[Phase],
    temperature: f64,
    pressure: f64,
    p0: f64,
    species_phase: &[usize],
    phase_stoichiometry: &[Vec<f64>],
) -> Result<Vec<f64>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();

    if gibbs.len() != m || element_totals.len() != elements.ncols() {
        return Err(ReactionExtentError::DimensionMismatch(
            "canonical residual input dimensions are inconsistent".to_string(),
        ));
    }
    validate_logmole_system_dimensions(
        reactions,
        elements,
        species_phase,
        phases.len(),
        phase_stoichiometry,
    )?;
    validate_residual_conditions(temperature, pressure, p0)?;
    if log_moles.len() != m {
        return Err(ReactionExtentError::ResidualEvaluation(
            "log-mole vector length mismatch".to_string(),
        ));
    }

    let moles = compute_species_moles(log_moles).map_err(|error| {
        ReactionExtentError::ResidualEvaluation(format!(
            "cannot reconstruct species moles for residual evaluation: {error}"
        ))
    })?;
    let mut phase_moles = vec![0.0; phases.len()];
    for (species, &phase) in species_phase.iter().enumerate() {
        phase_moles[phase] += moles[species];
    }
    for (index, &total) in phase_moles.iter().enumerate() {
        if !total.is_finite() || total <= 0.0 {
            return Err(ReactionExtentError::InvalidNPhase {
                index,
                value: total,
            });
        }
    }

    let phase_offsets = phase_activity_models(phases)
        .into_iter()
        .map(|model| model.log_phase_offset(pressure, p0))
        .collect::<Result<Vec<_>, _>>()?;
    let mut residual = vec![0.0; m];
    let rt = R * temperature;

    for reaction in 0..r {
        let mut sum_log_moles = 0.0;
        let mut sum_log_phase_moles = 0.0;
        let mut sum_phase_offsets = 0.0;
        let mut standard_gibbs = 0.0;
        for species in 0..m {
            let coefficient = reactions[(species, reaction)];
            if coefficient != 0.0 {
                sum_log_moles += coefficient * log_moles[species];
                let species_gibbs = gibbs[species](temperature);
                if !species_gibbs.is_finite() {
                    return Err(ReactionExtentError::InvalidDG0 {
                        species_index: species,
                        dg0: species_gibbs,
                        temperature,
                    });
                }
                standard_gibbs += coefficient * species_gibbs;
            }
        }
        for phase in 0..phases.len() {
            let coefficient = phase_stoichiometry[reaction][phase];
            if coefficient != 0.0 {
                sum_log_phase_moles += coefficient * phase_moles[phase].ln();
                sum_phase_offsets += coefficient * phase_offsets[phase];
            }
        }
        residual[reaction] =
            sum_log_moles - sum_log_phase_moles + sum_phase_offsets - (-standard_gibbs / rt);
        if !residual[reaction].is_finite() {
            return Err(ReactionExtentError::ResidualEvaluation(format!(
                "reaction {reaction} produced a non-finite residual"
            )));
        }
    }

    for element in 0..elements.ncols() {
        let reconstructed_total = (0..m)
            .map(|species| elements[(species, element)] * moles[species])
            .sum::<f64>();
        residual[r + element] = reconstructed_total - element_totals[element];
    }

    Ok(residual)
}

/// Deprecated threshold-based residual retained only for compatibility with
/// experiments outside the canonical equilibrium path.
#[deprecated(note = "use evaluate_equilibrium_logmole_residual instead")]
pub fn equilibrium_logmole_residual2(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs: Vec<GibbsFn>,      // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    subs_eps: f64,
    phase_eps: f64,
) -> Result<Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    let species_phase = species_to_phase_map(&phases, m)?;
    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    let rt = R * temperature;

    Ok(Box::new(move |y: &[f64]| {
        if y.len() != m {
            return Err(ReactionExtentError::ResidualEvaluation(
                "log-mole vector length mismatch".to_string(),
            ));
        }

        // --- species moles ---
        let mut n = vec![0.0; m];
        let mut active = vec![false; m];
        for i in 0..m {
            n[i] = y[i].exp();
            active[i] = species_active(n[i], subs_eps);
        }

        // --- phase totals ---
        let mut n_phase = vec![0.0; phases.len()];
        let mut phase_active_mask = vec![false; phases.len()];
        for i in 0..m {
            if active[i] {
                n_phase[species_phase[i]] += n[i];
            }
        }
        for j in 0..phases.len() {
            phase_active_mask[j] = phase_active(n_phase[j], phase_eps);
        }

        // --- phase offsets ---
        let phi = phase_activity_models(&phases)
            .into_iter()
            .map(|model| model.log_phase_offset(pressure, p0))
            .collect::<Result<Vec<_>, _>>()?;

        let mut f = vec![0.0; m];

        // -------------------------
        // Reaction equilibrium eqs
        // -------------------------
        for k in 0..r {
            let mut sum_ln_n = 0.0;
            let mut sum_ln_n_phase = 0.0;
            let mut sum_phi = 0.0;
            let mut dg0 = 0.0;

            for i in 0..m {
                if !active[i] {
                    continue;
                }
                let nu = reactions[(i, k)];
                if nu != 0.0 {
                    sum_ln_n += nu * y[i];
                    dg0 += nu * gibbs[i](temperature);
                }
            }

            for j in 0..phases.len() {
                if !phase_active_mask[j] {
                    continue;
                }
                let dnk = delta_n[k][j];
                if dnk != 0.0 {
                    sum_ln_n_phase += dnk * n_phase[j].ln();
                    sum_phi += dnk * phi[j];
                }
            }

            let ln_k = -dg0 / rt;
            f[k] = sum_ln_n - sum_ln_n_phase + sum_phi - ln_k;
        }

        // -------------------------
        // Element balance eqs
        // -------------------------
        for el in 0..e {
            let mut sum = 0.0;
            for i in 0..m {
                if active[i] {
                    sum += elements[(i, el)] * n[i];
                }
            }
            f[r + el] = sum - element_totals[el];
        }

        // -------------------------
        // Inactive species constraints
        // -------------------------
        for i in 0..m {
            if !active[i] {
                f[i] = y[i]; // push to -∞
            }
        }

        Ok(f)
    }))
}

/// Creates a scaled residual closure for improved numerical conditioning.
///
/// The scale vector must have one entry per residual row. Scaling is applied
/// after residual evaluation and before backend acceptance, so the physical
/// equations remain unchanged.
pub fn scaled_residual(
    f: Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>,
    scale: Vec<f64>,
) -> Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>> {
    Box::new(move |xi: &[f64]| scale_residual_rows(f(xi)?, &scale))
}

/// Applies a validated row scale without rebuilding a residual closure.
///
/// This helper is intentionally pure: it does not inspect solver state or
/// modify any cached problem data.
pub fn scale_residual_rows(
    mut residual: Vec<f64>,
    scale: &[f64],
) -> Result<Vec<f64>, ReactionExtentError> {
    if residual.len() != scale.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "residual has {} rows but scale has {} entries",
            residual.len(),
            scale.len(),
        )));
    }

    for (index, (value, factor)) in residual.iter_mut().zip(scale).enumerate() {
        if *factor <= 0.0 || !factor.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "residual_scale",
                message: format!("scale[{index}] must be finite and positive, got {factor}"),
            });
        }
        *value /= factor;
    }
    Ok(residual)
}

/// Evaluates the canonical analytical Jacobian for one log-moles state.
///
/// This is the pure counterpart to [`evaluate_equilibrium_logmole_residual`].
/// New backend adapters should use it instead of carrying the historical
/// threshold parameters required by [`equilibrium_logmole_jacobian`].
pub fn evaluate_equilibrium_logmole_jacobian(
    log_moles: &[f64],
    reactions: &DMatrix<f64>,
    elements: &DMatrix<f64>,
    species_phase: &[usize],
    phase_stoichiometry: &[Vec<f64>],
    phase_count: usize,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    equilibrium_logmole_jacobian(
        log_moles,
        reactions,
        elements,
        species_phase,
        phase_stoichiometry,
        phase_count,
        0.0,
        0.0,
    )
}

/// Computes the canonical Jacobian for the log-mole equilibrium residual.
///
/// The Jacobian follows the same row ordering as
/// [`evaluate_equilibrium_logmole_residual`]:
/// - reaction rows first,
/// - element-balance rows second.
///
/// Columns correspond to the log-mole coordinates in the prepared problem's
/// species ordering.
pub fn equilibrium_logmole_jacobian(
    y: &[f64],
    reactions: &DMatrix<f64>, // m × r
    elements: &DMatrix<f64>,  // m × E
    species_phase: &[usize],
    delta_n: &[Vec<f64>],
    n_phases: usize,

    _subs_eps: f64,
    _phase_eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    validate_logmole_system_dimensions(reactions, elements, species_phase, n_phases, delta_n)?;
    if y.len() != m {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "log-mole vector has {} entries, expected {m}",
            y.len()
        )));
    }

    let n = compute_species_moles(y).map_err(|error| {
        ReactionExtentError::JacobianEvaluation(format!(
            "cannot reconstruct species moles for Jacobian evaluation: {error}"
        ))
    })?;

    let mut n_phase = vec![0.0; n_phases];
    for i in 0..m {
        n_phase[species_phase[i]] += n[i];
    }
    for (index, &total) in n_phase.iter().enumerate() {
        if !total.is_finite() || total <= 0.0 {
            return Err(ReactionExtentError::InvalidNPhase {
                index,
                value: total,
            });
        }
    }

    let mut jmat = DMatrix::<f64>::zeros(m, m);

    // --- equilibrium equations ---
    for k in 0..r {
        for i in 0..m {
            let nu = reactions[(i, k)];
            // A species with zero direct stoichiometric coefficient still
            // differentiates ln(N_phase) when it shares a phase with a
            // reacting species. Omitting this term breaks mixed-phase
            // Jacobians while leaving many all-reacting gas fixtures green.
            let phase = species_phase[i];
            let term = delta_n[k][phase] * n[i] / n_phase[phase];
            jmat[(k, i)] = nu - term;
        }
    }

    // --- element balances ---
    for el in 0..e {
        for i in 0..m {
            jmat[(r + el, i)] = elements[(i, el)] * n[i];
        }
    }

    Ok(jmat)
}

/// Deprecated threshold-based Jacobian paired with
/// [`equilibrium_logmole_residual2`].
#[deprecated(note = "use evaluate_equilibrium_logmole_jacobian instead")]
pub fn equilibrium_logmole_jacobian2(
    y: &[f64],
    reactions: &DMatrix<f64>, // m × r
    elements: &DMatrix<f64>,  // m × E
    species_phase: &[usize],
    delta_n: &[Vec<f64>],
    n_phases: usize,
    subs_eps: f64,
    phase_eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    let mut n = vec![0.0; m];
    let mut active = vec![false; m];
    for i in 0..m {
        n[i] = y[i].exp();
        active[i] = species_active(n[i], subs_eps);
    }

    let mut n_phase = vec![0.0; n_phases];
    let mut phase_active_mask = vec![false; n_phases];
    for i in 0..m {
        if active[i] {
            n_phase[species_phase[i]] += n[i];
        }
    }
    for j in 0..n_phases {
        phase_active_mask[j] = phase_active(n_phase[j], phase_eps);
    }

    let mut jmat = DMatrix::<f64>::zeros(m, m);

    // --- reaction equilibrium rows ---
    for k in 0..r {
        for i in 0..m {
            if !active[i] {
                continue;
            }
            let nu = reactions[(i, k)];

            let phase = species_phase[i];
            if phase_active_mask[phase] {
                let term = delta_n[k][phase] * n[i] / n_phase[phase];
                jmat[(k, i)] = nu - term;
            } else {
                jmat[(k, i)] = nu;
            }
        }
    }

    // --- element balances ---
    for el in 0..e {
        for i in 0..m {
            if active[i] {
                jmat[(r + el, i)] = elements[(i, el)] * n[i];
            }
        }
    }

    // --- inactive species diagonal ---
    for i in 0..m {
        if !active[i] {
            jmat[(i, i)] = 1.0;
        }
    }

    Ok(jmat)
}

/// Creates a scaled Jacobian closure for improved conditioning.
///
/// The same row scale used for the residual must be reused here so that the
/// backend sees a consistent linearization of the same equations.
pub fn scaled_jacobian(
    j: Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>>,
    scale: Vec<f64>,
) -> Box<dyn Fn(&[f64]) -> Result<DMatrix<f64>, ReactionExtentError>> {
    Box::new(move |xi: &[f64]| scale_jacobian_rows(j(xi)?, &scale))
}

/// Applies a validated row scale to an analytical Jacobian.
pub fn scale_jacobian_rows(
    mut jacobian: DMatrix<f64>,
    scale: &[f64],
) -> Result<DMatrix<f64>, ReactionExtentError> {
    if jacobian.nrows() != scale.len() {
        return Err(ReactionExtentError::DimensionMismatch(format!(
            "Jacobian has {} rows but scale has {} entries",
            jacobian.nrows(),
            scale.len(),
        )));
    }

    for k in 0..jacobian.nrows() {
        let s = scale[k];
        if s <= 0.0 || !s.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "residual_scale",
                message: format!("scale[{k}] must be finite and positive, got {s}"),
            });
        }
        for jcol in 0..jacobian.ncols() {
            jacobian[(k, jcol)] /= s;
        }
    }
    Ok(jacobian)
}
