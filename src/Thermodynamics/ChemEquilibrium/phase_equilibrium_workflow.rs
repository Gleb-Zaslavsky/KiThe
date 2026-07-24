//! Narrow public fixed-`P,T` phase-equilibrium workflow.
//!
//! This facade owns orchestration only: it joins validated resolved data,
//! physical inventory, numerical settings, bridge construction, and immutable
//! result publication. It does not duplicate residual construction or expose
//! the historical mutable solver as an alternative public engine.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::EquilibriumSolverSettings;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::MultiphaseInitialComposition;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, TraceSpeciesSeedPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_problem::{
    PhaseEquilibriumBuildRequest, SupportedPhaseModelPolicy, build_phase_equilibrium_problem,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;

/// Current solve mode exposed by the typed public facade.
///
/// Bounded phase-control orchestration remains an internal migration target.
/// The facade names the implemented fixed-declared-phase contract explicitly
/// instead of pretending that legacy mutable phase control is equivalent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PhaseEquilibriumSolveMode {
    /// Solve all declared phases as one immutable fixed active set.
    #[default]
    FixedDeclaredPhases,
}

/// Complete public input for one resolved fixed-pressure, fixed-temperature
/// equilibrium solve.
#[derive(Clone)]
pub struct ResolvedPhaseEquilibriumRequest<'a> {
    resolved: &'a ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    initial_composition: MultiphaseInitialComposition,
    trace_seed_policy: TraceSpeciesSeedPolicy,
    model_policy: SupportedPhaseModelPolicy,
    solver_settings: EquilibriumSolverSettings,
    solve_mode: PhaseEquilibriumSolveMode,
}

impl<'a> ResolvedPhaseEquilibriumRequest<'a> {
    /// Creates one explicit solver request. The resolved system remains
    /// borrowed and immutable for the full build/solve transaction.
    pub fn new(
        resolved: &'a ResolvedPhaseSystem,
        conditions: EquilibriumConditions,
        initial_composition: MultiphaseInitialComposition,
    ) -> Self {
        Self {
            resolved,
            conditions,
            initial_composition,
            trace_seed_policy: TraceSpeciesSeedPolicy::Absolute {
                floor: crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::DEFAULT_TRACE_MOLE_FLOOR,
            },
            model_policy: SupportedPhaseModelPolicy::default(),
            solver_settings: EquilibriumSolverSettings::default(),
            solve_mode: PhaseEquilibriumSolveMode::default(),
        }
    }

    /// Replaces only the numerical trace-coordinate policy.
    pub fn with_trace_seed_policy(mut self, policy: TraceSpeciesSeedPolicy) -> Self {
        self.trace_seed_policy = policy;
        self
    }

    /// Replaces only the supported physical-model policy.
    pub fn with_model_policy(mut self, policy: SupportedPhaseModelPolicy) -> Self {
        self.model_policy = policy;
        self
    }

    /// Replaces only numerical backend and acceptance settings.
    pub fn with_solver_settings(mut self, settings: EquilibriumSolverSettings) -> Self {
        self.solver_settings = settings;
        self
    }

    /// Selects the explicit orchestration mode.
    pub fn with_solve_mode(mut self, mode: PhaseEquilibriumSolveMode) -> Self {
        self.solve_mode = mode;
        self
    }

    /// Borrow the immutable resolved system.
    pub fn resolved(&self) -> &'a ResolvedPhaseSystem {
        self.resolved
    }

    /// Fixed thermodynamic conditions.
    pub fn conditions(&self) -> EquilibriumConditions {
        self.conditions
    }

    /// Physical initial inventory in canonical layout order.
    pub fn initial_composition(&self) -> &MultiphaseInitialComposition {
        &self.initial_composition
    }

    /// Numerical policy to be applied after the physical inventory is validated.
    pub fn solver_settings(&self) -> &EquilibriumSolverSettings {
        &self.solver_settings
    }
}

/// Builds and solves one resolved phase system through the canonical bridge.
///
/// All preparation stays transactional: a failed lookup, validation, or
/// backend attempt returns an error without modifying the supplied resolved
/// phase system or publishing a partial result.
pub fn solve_resolved_pt(
    request: ResolvedPhaseEquilibriumRequest<'_>,
) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
    match request.solve_mode {
        PhaseEquilibriumSolveMode::FixedDeclaredPhases => {
            let bundle = build_phase_equilibrium_problem(PhaseEquilibriumBuildRequest::new(
                request.resolved,
                request.conditions,
                request.initial_composition,
                request.trace_seed_policy,
                request.model_policy,
            )?)?;
            bundle
                .solve_with(|settings| *settings = request.solver_settings)
                .and_then(|bundle| bundle.into_multiphase_solution())
        }
    }
}
