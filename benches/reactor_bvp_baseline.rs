//! PERF-BVP-01: matched symbolic BVP routes at increasing system sizes.
//!
//! Criterion measures cold problem construction through accepted publication
//! with telemetry Off. A separate Detailed pass reports fixture parsing,
//! problem validation, solver preparation, native solve, and route telemetry.
//! The analytic profiles are dimensionless solver controls, not combustion
//! validation fixtures; physical reactor acceptance lives in BVP story tests.
//!
//! Run with:
//! cargo bench --bench reactor_bvp_baseline --no-default-features -j 2 -- --noplot
//! KITHE_RUN_BVP_AOT_BENCH=1 enables the optional TCC AOT/Lambdify axis.

use std::collections::HashMap;
use std::hint::black_box;
use std::path::Path;
use std::time::{Duration, Instant};

use criterion::{BatchSize, BenchmarkId, Criterion, criterion_group, criterion_main};
use nalgebra::DMatrix;
use tabled::settings::Style;
use tabled::{Table, Tabled};

use KiThe::SimpleReactorsBVP::SimpleReactorBVP::{
    FastElemReact, SimpleReactorTask, SymbolicRhsAssemblyPolicy,
};
use KiThe::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use KiThe::SimpleReactorsBVP::solver_backend::{
    ReactorBvpAotCompiler, ReactorBvpMatrixBackend, ReactorBvpSolverConfig,
    ReactorBvpSymbolicBackend,
};
use RustedSciThe::numerical::BVP_CorrectionNewton::{
    CorrectionNewtonMatrixLayout, CorrectionNewtonOptions, CorrectionNewtonSymbolicExecution,
    CorrectionNewtonSymbolicFrontend, CorrectionNewtonTelemetryMode,
};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{AdaptiveGridConfig, SolverParams};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_frozen::FrozenSolverOptions;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;
use RustedSciThe::numerical::BVP_Damp::telemetry::BvpTelemetryMode;
use RustedSciThe::numerical::BVP_api::{
    BvpBoundaryCondition, BvpBoundarySide, BvpProblem, BvpSciSparseFactorization, BvpSolver,
    BvpSolverConfig, BvpSolverTelemetrySnapshot,
};
use RustedSciThe::numerical::BVP_sci::new::{
    BvpSciAssembly, BvpSciExecution, BvpSciMatrixLayout, BvpSciOptions, BvpSciTelemetryMode,
};
use RustedSciThe::symbolic::bvp::telemetry::BvpLambdifyExecutionPolicy;
use RustedSciThe::symbolic::codegen::codegen_provider_api::MatrixBackend;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend;
use RustedSciThe::symbolic::symbolic_ivp_generated::SymbolicIvpGeneratedBackendConfig;

mod large_mechanism_fixture {
    use KiThe as root;
    include!("support/large_mechanism_fixture.rs");
}

const DIAGNOSTIC_REPETITIONS: usize = 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Family {
    CorrectionNewton,
    Damped,
    Frozen,
    Naive,
    Sci,
}

impl Family {
    fn label(self) -> &'static str {
        match self {
            Self::CorrectionNewton => "CorrectionNewton",
            Self::Damped => "Damped",
            Self::Frozen => "Frozen",
            Self::Naive => "Naive",
            Self::Sci => "BVP_sci",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Frontend {
    ExprLegacy,
    AtomView,
}

impl Frontend {
    fn label(self) -> &'static str {
        match self {
            Self::ExprLegacy => "ExprLegacy",
            Self::AtomView => "AtomView",
        }
    }

    fn damped(self) -> ReactorBvpSymbolicBackend {
        match self {
            Self::ExprLegacy => ReactorBvpSymbolicBackend::ExprLegacy,
            Self::AtomView => ReactorBvpSymbolicBackend::AtomView,
        }
    }

    fn sci(self) -> BvpSciAssembly {
        match self {
            Self::ExprLegacy => BvpSciAssembly::ExprLegacy,
            Self::AtomView => BvpSciAssembly::AtomViewNative,
        }
    }

    fn frozen(self) -> BvpSymbolicAssemblyBackend {
        match self {
            Self::ExprLegacy => BvpSymbolicAssemblyBackend::ExprLegacy,
            Self::AtomView => BvpSymbolicAssemblyBackend::AtomView,
        }
    }

    fn correction_newton(self) -> CorrectionNewtonSymbolicFrontend {
        match self {
            Self::ExprLegacy => CorrectionNewtonSymbolicFrontend::ExprLegacy,
            Self::AtomView => CorrectionNewtonSymbolicFrontend::AtomView,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Layout {
    Banded,
    Sparse,
    Dense,
}

impl Layout {
    fn label(self, family: Family) -> &'static str {
        match (family, self) {
            (Family::CorrectionNewton, Self::Banded) => "BandedAuto",
            (Family::CorrectionNewton, Self::Sparse) => "Sparse",
            (Family::CorrectionNewton, Self::Dense) => "Dense",
            (Family::Damped | Family::Frozen | Family::Naive, Self::Banded) => "Banded",
            (Family::Damped | Family::Frozen | Family::Naive, Self::Sparse) => "Sparse",
            (Family::Frozen | Family::Naive, Self::Dense) => "Dense",
            (Family::Sci, Self::Banded) => "BandedAuto",
            (Family::Sci, Self::Sparse) => "Sparse",
            (Family::Sci, Self::Dense) => "Dense",
            (Family::Damped, Self::Dense) => unreachable!("Dense is not a Damped route"),
        }
    }

    fn sci(self) -> BvpSciMatrixLayout {
        match self {
            Self::Banded => BvpSciMatrixLayout::BandedAuto,
            Self::Sparse => BvpSciMatrixLayout::Sparse,
            Self::Dense => BvpSciMatrixLayout::Dense,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Scenario {
    name: &'static str,
    states: usize,
    intervals: usize,
}

const SCENARIOS: [Scenario; 3] = [
    Scenario {
        name: "scalar-small",
        states: 1,
        intervals: 24,
    },
    Scenario {
        name: "eight-state-medium",
        states: 8,
        intervals: 64,
    },
    Scenario {
        name: "twenty-four-state-large",
        states: 24,
        intervals: 128,
    },
];

#[derive(Clone, Copy, Debug)]
struct Route {
    family: Family,
    frontend: Frontend,
    layout: Layout,
    sci_factorization: BvpSciSparseFactorization,
}

impl Route {
    fn new(family: Family, frontend: Frontend, layout: Layout) -> Self {
        Self {
            family,
            frontend,
            layout,
            sci_factorization: BvpSciSparseFactorization::StandardFaer,
        }
    }

    fn label(self) -> String {
        let label = format!(
            "{}/{}/{}",
            self.family.label(),
            self.frontend.label(),
            self.layout.label(self.family)
        );
        match self.sci_factorization {
            BvpSciSparseFactorization::StandardFaer => label,
            BvpSciSparseFactorization::SuperLu => format!("{label}/SuperLU"),
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum CallbackPolicy {
    Sequential,
    Auto,
    ForcedParallel,
}

impl CallbackPolicy {
    fn label(self) -> &'static str {
        match self {
            Self::Sequential => "Sequential",
            Self::Auto => "Auto",
            Self::ForcedParallel => "ForcedParallel",
        }
    }

    fn to_rst(self) -> BvpLambdifyExecutionPolicy {
        match self {
            Self::Sequential => BvpLambdifyExecutionPolicy::Sequential,
            Self::Auto => BvpLambdifyExecutionPolicy::Auto { min_work: 0 },
            Self::ForcedParallel => BvpLambdifyExecutionPolicy::Parallel { min_work: 0 },
        }
    }
}

#[derive(Clone, Debug)]
struct Sample {
    equation_build: Duration,
    problem_validation: Duration,
    solver_preparation: Duration,
    native_solve: Duration,
    total: Duration,
    mesh_nodes: usize,
    max_analytic_error: f64,
    values: Vec<f64>,
    telemetry: BvpSolverTelemetrySnapshot,
}

#[derive(Tabled)]
struct MacroRow {
    scenario: String,
    route: String,
    execution: String,
    states: usize,
    nodes: usize,
    max_error: String,
    preparation_ms: String,
    solve_ms: String,
    total_ms: String,
}

#[derive(Tabled)]
struct StageRow {
    scenario: String,
    route: String,
    execution: String,
    stage: String,
    median_ms: String,
}

#[derive(Tabled)]
struct NativeRow {
    scenario: String,
    route: String,
    execution: String,
    residual_or_rhs_calls: String,
    jacobian_calls_or_builds: String,
    linear_assemblies: String,
    factorizations: String,
    linear_solves: String,
    collocation_evaluations: String,
    mesh_refinements: String,
    residual_or_rhs_ms: String,
    jacobian_ms: String,
    linear_assembly_ms: String,
    factorization_ms: String,
    linear_solve_ms: String,
    collocation_ms: String,
    mesh_refinement_ms: String,
    solver_internal_ms: String,
    full_solve_ms: String,
    callback_parallel_dispatches: String,
    callback_sequential_dispatches: String,
    superlu_factor_calls: String,
    superlu_factor_ms: String,
    superlu_solve_ms: String,
}

#[derive(Tabled)]
struct PreparationStageRow {
    fixture: String,
    assembly_policy: String,
    stage: String,
    median_ms: String,
}

#[derive(Tabled)]
struct CallbackPolicyRow {
    scenario: String,
    solver: String,
    policy: String,
    max_profile_delta: String,
    parallel_dispatches: u64,
    sequential_dispatches: u64,
    callback_ms: String,
    median_solve_ms: String,
}

#[derive(Tabled)]
struct ComparisonRow {
    scenario: String,
    route: String,
    lambdify_total_ms: String,
    aot_total_ms: String,
    aot_delta_pct: String,
}

fn routes() -> Vec<Route> {
    let mut routes = Vec::new();
    for frontend in [Frontend::ExprLegacy, Frontend::AtomView] {
        for layout in [Layout::Dense, Layout::Sparse, Layout::Banded] {
            routes.push(Route::new(Family::CorrectionNewton, frontend, layout));
        }
        for layout in [Layout::Banded, Layout::Sparse] {
            routes.push(Route::new(Family::Damped, frontend, layout));
        }
        for layout in [Layout::Dense, Layout::Sparse, Layout::Banded] {
            routes.push(Route::new(Family::Sci, frontend, layout));
        }
        for family in [Family::Frozen, Family::Naive] {
            for layout in [Layout::Dense, Layout::Sparse, Layout::Banded] {
                routes.push(Route::new(family, frontend, layout));
            }
        }
    }
    if superlu_bridge_path().is_some() {
        for frontend in [Frontend::ExprLegacy, Frontend::AtomView] {
            let mut route = Route::new(Family::Sci, frontend, Layout::Sparse);
            route.sci_factorization = BvpSciSparseFactorization::SuperLu;
            routes.push(route);
        }
    }
    routes
}

fn superlu_bridge_path() -> Option<std::path::PathBuf> {
    ["SUPERLU_LIBRARY", "SUPERLU_EXPERIMENTAL_LIBRARY"]
        .into_iter()
        .filter_map(|name| std::env::var_os(name).map(std::path::PathBuf::from))
        .find(|path| path.is_file())
}

fn large_reactor_setup_fixture(policy: SymbolicRhsAssemblyPolicy) -> SimpleReactorTask {
    let species = (1..=40).map(|n| format!("C{n}")).collect::<Vec<_>>();
    let reactions = (1..40)
        .map(|n| FastElemReact {
            eq: if n == 1 {
                "C1+C1=>C2".to_string()
            } else {
                format!("C{n}+C1=>C{}", n + 1)
            },
            A: 1.0e3,
            n: 0.0,
            E: 1.0e4,
            Q: 0.0,
        })
        .collect::<Vec<_>>();
    let mut reactor = SimpleReactorTask::new().with_symbolic_rhs_assembly_policy(policy);
    reactor
        .fast_react_set(reactions)
        .expect("large setup fixture reactions should parse");
    let boundary = species
        .iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), if i == 0 { 1.0 } else { 0.0 }))
        .chain([("T".to_string(), 700.0)])
        .collect::<HashMap<_, _>>();
    let diffusion = species
        .iter()
        .map(|name| (name.clone(), 1.0e-5))
        .collect::<HashMap<_, _>>();
    reactor.set_parameters(
        vec![0.0; 39],
        101_325.0,
        700.0,
        1_000.0,
        boundary,
        0.1,
        diffusion,
        0.01,
        ScalingConfig::new(600.0, 1.0e-3, 600.0),
    );
    reactor.M = 0.012;
    reactor
}

fn reactor_setup_policy_rows() -> Vec<PreparationStageRow> {
    let mut rows = Vec::new();
    let mut reference: Option<Vec<String>> = None;
    for policy in [
        SymbolicRhsAssemblyPolicy::Sequential,
        SymbolicRhsAssemblyPolicy::Auto,
        SymbolicRhsAssemblyPolicy::Parallel,
    ] {
        let mut elapsed = Vec::new();
        let mut equations = Vec::new();
        for _ in 0..DIAGNOSTIC_REPETITIONS {
            let mut reactor = large_reactor_setup_fixture(policy);
            let started = Instant::now();
            reactor
                .setup_bvp()
                .expect("40-species KiThe setup should succeed before RST handoff");
            elapsed.push(ms(started.elapsed()));
            equations = reactor
                .solver
                .eq_system
                .iter()
                .map(ToString::to_string)
                .collect();
        }
        if let Some(reference) = &reference {
            assert_eq!(
                reference, &equations,
                "RHS assembly policy changed equations"
            );
        } else {
            reference = Some(equations);
        }
        rows.push(PreparationStageRow {
            fixture: "40 species / 39 reactions".into(),
            assembly_policy: format!("{policy:?}"),
            stage: "KiThe setup before RST handoff".into(),
            median_ms: format!("{:.4}", median(&mut elapsed)),
        });
    }
    rows
}

fn is_aot_enabled() -> bool {
    matches!(
        std::env::var("KITHE_RUN_BVP_AOT_BENCH").ok().as_deref(),
        Some("1") | Some("true") | Some("TRUE") | Some("yes") | Some("YES")
    ) && std::process::Command::new("tcc")
        .arg("-v")
        .output()
        .is_ok_and(|output| output.status.success())
}

fn supports_aot(route: Route) -> bool {
    if route.sci_factorization == BvpSciSparseFactorization::SuperLu {
        return false;
    }
    match route.family {
        Family::CorrectionNewton | Family::Damped | Family::Sci => true,
        Family::Frozen | Family::Naive => {
            route.frontend == Frontend::AtomView
                && matches!(route.layout, Layout::Sparse | Layout::Banded)
        }
    }
}

fn equations(scenario: Scenario) -> (Vec<Expr>, Vec<String>, Vec<f64>) {
    let rates = (0..scenario.states)
        .map(|index| 0.25 + 0.035 * index as f64)
        .collect::<Vec<_>>();
    let names = (0..scenario.states)
        .map(|index| format!("y{index}"))
        .collect::<Vec<_>>();
    let equations = rates
        .iter()
        .zip(&names)
        .map(|(rate, name)| Expr::parse_expression(&format!("-({rate:.17}) * {name}")))
        .collect();
    (equations, names, rates)
}

fn build_route_config(
    route: Route,
    initial_guess: DMatrix<f64>,
    telemetry: bool,
    aot_output: Option<&Path>,
    callback_policy: Option<CallbackPolicy>,
) -> BvpSolverConfig {
    match route.family {
        Family::CorrectionNewton => {
            let matrix_layout = match route.layout {
                Layout::Dense => CorrectionNewtonMatrixLayout::Dense,
                Layout::Sparse => CorrectionNewtonMatrixLayout::Sparse,
                Layout::Banded => CorrectionNewtonMatrixLayout::BandedAuto,
            };
            let mut options = CorrectionNewtonOptions {
                symbolic_frontend: route.frontend.correction_newton(),
                symbolic_execution: if aot_output.is_some() {
                    CorrectionNewtonSymbolicExecution::Aot
                } else {
                    CorrectionNewtonSymbolicExecution::Lambdify
                },
                matrix_layout,
                absolute_tolerance: 1e-8,
                default_relative_tolerance: 1e-8,
                max_iterations: 100,
                telemetry_mode: if telemetry {
                    CorrectionNewtonTelemetryMode::Timings
                } else {
                    CorrectionNewtonTelemetryMode::Off
                },
                ..CorrectionNewtonOptions::default()
            };
            options.bounds = bounds_for_damped(initial_guess.nrows());
            options.relative_tolerances = rel_tolerance_for_damped(initial_guess.nrows());
            if let Some(path) = aot_output {
                options = options.with_aot_config(
                    SymbolicIvpGeneratedBackendConfig::build_if_missing_release(path).with_c_tcc(),
                );
            }
            BvpSolverConfig::CorrectionNewton {
                options,
                initial_guess,
            }
        }
        Family::Damped => {
            let matrix = match route.layout {
                Layout::Banded => ReactorBvpMatrixBackend::Banded,
                Layout::Sparse => ReactorBvpMatrixBackend::Sparse,
                Layout::Dense => unreachable!("Dense is not a Damped capability"),
            };
            let config = match (route.layout, aot_output.is_some()) {
                (Layout::Banded, true) => {
                    ReactorBvpSolverConfig::banded_aot(ReactorBvpAotCompiler::CTcc)
                }
                (Layout::Sparse, true) => {
                    ReactorBvpSolverConfig::sparse_aot(ReactorBvpAotCompiler::CTcc)
                }
                (_, false) => ReactorBvpSolverConfig::default_lambdify(),
                (Layout::Dense, _) => unreachable!("Dense is not a Damped capability"),
            }
            .with_matrix_backend(matrix)
            .with_symbolic_backend(route.frontend.damped());
            let telemetry = if telemetry {
                BvpTelemetryMode::Detailed
            } else {
                BvpTelemetryMode::Off
            };
            let mut options = config
                .to_rusted_options()
                .expect("Damped backend configuration should be supported")
                .with_abs_tolerance(1e-8)
                .with_max_iterations(100)
                .with_rel_tolerance(rel_tolerance_for_damped(initial_guess.nrows()))
                .with_bounds(bounds_for_damped(initial_guess.nrows()))
                .with_bvp_telemetry_mode(telemetry);
            if let Some(policy) = callback_policy {
                let mut generated = options.generated_backend_config.clone();
                generated.lambdify_execution_policy = policy.to_rst();
                generated.lambdify_telemetry_mode =
                    if matches!(telemetry, BvpTelemetryMode::Detailed) {
                        RustedSciThe::symbolic::bvp::telemetry::BvpLambdifyTelemetryMode::Detailed
                    } else {
                        RustedSciThe::symbolic::bvp::telemetry::BvpLambdifyTelemetryMode::Off
                    };
                options = options.with_generated_backend_config(generated);
            }
            BvpSolverConfig::Damped {
                options,
                initial_guess,
            }
        }
        Family::Sci => {
            let execution = if aot_output.is_some() {
                BvpSciExecution::Aot
            } else {
                BvpSciExecution::Lambdify
            };
            let options = BvpSciOptions::default()
                .with_execution(execution)
                .with_assembly(route.frontend.sci())
                .with_matrix_layout(route.layout.sci())
                .with_tolerance(1e-6)
                .with_boundary_tolerance(1e-7)
                .with_limits(4 * initial_guess.ncols(), 20, 6)
                .with_execution_policy(
                    callback_policy.map_or(BvpLambdifyExecutionPolicy::Sequential, |policy| {
                        policy.to_rst()
                    }),
                )
                .with_telemetry(if telemetry {
                    BvpSciTelemetryMode::Detailed
                } else {
                    BvpSciTelemetryMode::Off
                });
            let aot_config = aot_output.map(|path| {
                SymbolicIvpGeneratedBackendConfig::build_if_missing_release(path).with_c_tcc()
            });
            BvpSolverConfig::Sci {
                options,
                sparse_factorization: route.sci_factorization,
                aot_config,
                initial_guess,
            }
        }
        Family::Frozen | Family::Naive => {
            let matrix = match route.layout {
                Layout::Dense => MatrixBackend::Dense,
                Layout::Sparse => MatrixBackend::SparseCol,
                Layout::Banded => MatrixBackend::Banded,
            };
            let mut options = match (route.layout, aot_output.is_some()) {
                (Layout::Dense, false) => FrozenSolverOptions::dense_frozen(),
                (Layout::Sparse, false) => FrozenSolverOptions::sparse_frozen(),
                (Layout::Banded, false) => {
                    FrozenSolverOptions::banded_frozen().with_banded_lambdify()
                }
                (Layout::Sparse, true) => {
                    FrozenSolverOptions::sparse_frozen().with_sparse_atomview_c_tcc()
                }
                (Layout::Banded, true) => {
                    FrozenSolverOptions::banded_frozen().with_banded_atomview_c_tcc()
                }
                (Layout::Dense, true) => unreachable!("Dense has no Frozen AOT route"),
            }
            .with_tolerance(1e-8)
            .with_max_iterations(500)
            .with_bvp_telemetry_mode(if telemetry {
                BvpTelemetryMode::Detailed
            } else {
                BvpTelemetryMode::Off
            })
            .with_symbolic_assembly_backend(route.frontend.frozen())
            .with_matrix_backend(matrix);
            if route.family == Family::Naive {
                options.strategy = "Naive".to_string();
                options.strategy_params = None;
            }
            BvpSolverConfig::Frozen {
                options,
                initial_guess,
            }
        }
    }
}

fn rel_tolerance_for_damped(states: usize) -> HashMap<String, f64> {
    (0..states)
        .map(|index| (format!("y{index}"), 1e-8))
        .collect()
}

fn bounds_for_damped(states: usize) -> HashMap<String, (f64, f64)> {
    (0..states)
        .map(|index| (format!("y{index}"), (0.0, 1.1)))
        .collect()
}

fn run_once(
    scenario: Scenario,
    route: Route,
    telemetry: bool,
    aot_output: Option<&Path>,
) -> Sample {
    run_once_with_policy(scenario, route, telemetry, aot_output, None)
}

fn run_once_with_policy(
    scenario: Scenario,
    route: Route,
    telemetry: bool,
    aot_output: Option<&Path>,
    callback_policy: Option<CallbackPolicy>,
) -> Sample {
    let total_started = Instant::now();
    let started = Instant::now();
    let (equations, names, rates) = equations(scenario);
    let equation_build = started.elapsed();

    let mesh = (0..=scenario.intervals)
        .map(|index| index as f64 / scenario.intervals as f64)
        .collect::<Vec<_>>();
    let conditions = (0..scenario.states)
        .map(|state| BvpBoundaryCondition {
            state,
            side: BvpBoundarySide::Left,
            value: 1.0,
        })
        .collect();
    let started = Instant::now();
    let problem = BvpProblem::new(equations, names, "x", mesh.clone(), conditions)
        .expect("analytic fixture should form a valid BVP");
    let problem_validation = started.elapsed();

    let guess_columns = match route.family {
        Family::CorrectionNewton | Family::Damped => scenario.intervals,
        Family::Frozen | Family::Naive => scenario.intervals,
        Family::Sci => scenario.intervals + 1,
    };
    let initial_guess = DMatrix::from_fn(scenario.states, guess_columns, |state, column| {
        let node = if route.family == Family::CorrectionNewton {
            column + 1
        } else {
            column
        };
        let x = node as f64 / scenario.intervals as f64;
        (-rates[state] * x).exp()
    });
    let config = build_route_config(route, initial_guess, telemetry, aot_output, callback_policy);
    let started = Instant::now();
    let mut solver = BvpSolver::new(problem, config).unwrap_or_else(|error| {
        panic!(
            "{} / {} preparation failed: {error}",
            scenario.name,
            route.label()
        )
    });
    let solver_preparation = started.elapsed();
    let started = Instant::now();
    solver.solve().unwrap_or_else(|error| {
        panic!(
            "{} / {} solve failed: {error}",
            scenario.name,
            route.label()
        )
    });
    let native_solve = started.elapsed();
    let view = solver
        .solution_view()
        .expect("successful solve retains output");
    let mesh_points = view.mesh().len();
    let mut values = Vec::with_capacity(mesh_points * scenario.states);
    let mut max_analytic_error: f64 = 0.0;
    for node in 0..mesh_points {
        for (state, rate) in rates.iter().enumerate() {
            let value = view
                .value(node, state)
                .expect("solution view index is valid");
            assert!(value.is_finite(), "solution contains non-finite values");
            assert!(
                (0.0..=1.0 + 1e-10).contains(&value),
                "dimensionless state escaped [0,1]: {value}"
            );
            let x = view.mesh()[node];
            max_analytic_error = max_analytic_error.max((value - (-rate * x).exp()).abs());
            values.push(value);
        }
    }
    assert!(
        max_analytic_error <= 0.05,
        "{} / {} analytic error {max_analytic_error:e} exceeds smoke envelope",
        scenario.name,
        route.label()
    );
    let telemetry_snapshot = solver.telemetry_snapshot();
    if telemetry && route.sci_factorization == BvpSciSparseFactorization::SuperLu {
        let BvpSolverTelemetrySnapshot::Sci(snapshot) = &telemetry_snapshot else {
            unreachable!("SuperLU is only selected for BVP_sci")
        };
        assert!(
            snapshot.superlu_fallbacks == 0,
            "SuperLU unexpectedly fell back to Faer"
        );
        assert!(
            snapshot.superlu.native_factorization_calls > 0,
            "SuperLU was requested but no native SuperLU factorization was observed"
        );
    }
    Sample {
        equation_build,
        problem_validation,
        solver_preparation,
        native_solve,
        total: total_started.elapsed(),
        mesh_nodes: mesh_points,
        max_analytic_error,
        values,
        telemetry: telemetry_snapshot,
    }
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn ms(value: Duration) -> f64 {
    value.as_secs_f64() * 1000.0
}

fn optional_ms(value: Option<f64>) -> String {
    value.map_or_else(
        || "n/a".into(),
        |milliseconds| format!("{:.4}", milliseconds.max(0.0)),
    )
}

fn collect_native_rows(
    scenario: Scenario,
    route: Route,
    execution: &str,
    telemetry: &BvpSolverTelemetrySnapshot,
) -> Vec<NativeRow> {
    let mut rows = Vec::new();
    match telemetry {
        BvpSolverTelemetrySnapshot::CorrectionNewton(Some(snapshot)) => {
            rows.push(NativeRow {
                scenario: scenario.name.into(),
                route: route.label(),
                execution: execution.into(),
                residual_or_rhs_calls: snapshot.residual_evaluations.to_string(),
                jacobian_calls_or_builds: snapshot.jacobian_evaluations.to_string(),
                linear_assemblies: "n/a".into(),
                factorizations: snapshot.factorizations.to_string(),
                linear_solves: snapshot.linear_solves.to_string(),
                collocation_evaluations: "n/a".into(),
                mesh_refinements: "n/a; fixed mesh".into(),
                residual_or_rhs_ms: format!("{:.4}", ms(snapshot.residual_callback_time)),
                jacobian_ms: format!("{:.4}", ms(snapshot.jacobian_callback_time)),
                linear_assembly_ms: format!("{:.4}", ms(snapshot.jacobian_assembly_time)),
                factorization_ms: format!("{:.4}", ms(snapshot.factorization_time)),
                linear_solve_ms: format!("{:.4}", ms(snapshot.linear_solve_time)),
                collocation_ms: "n/a".into(),
                mesh_refinement_ms: "n/a; fixed mesh".into(),
                solver_internal_ms: "n/a; no aggregate timer".into(),
                full_solve_ms: "n/a; use outer solve wall".into(),
                callback_parallel_dispatches: "n/a".into(),
                callback_sequential_dispatches: "n/a".into(),
                superlu_factor_calls: "n/a".into(),
                superlu_factor_ms: "n/a".into(),
                superlu_solve_ms: "n/a".into(),
            });
        }
        BvpSolverTelemetrySnapshot::CorrectionNewton(None) => {
            rows.push(NativeRow {
                scenario: scenario.name.into(),
                route: route.label(),
                execution: execution.into(),
                residual_or_rhs_calls: "n/a; telemetry off".into(),
                jacobian_calls_or_builds: "n/a; telemetry off".into(),
                linear_assemblies: "n/a".into(),
                factorizations: "n/a; telemetry off".into(),
                linear_solves: "n/a; telemetry off".into(),
                collocation_evaluations: "n/a".into(),
                mesh_refinements: "n/a; fixed mesh".into(),
                residual_or_rhs_ms: "n/a; telemetry off".into(),
                jacobian_ms: "n/a; telemetry off".into(),
                linear_assembly_ms: "n/a; telemetry off".into(),
                factorization_ms: "n/a; telemetry off".into(),
                linear_solve_ms: "n/a; telemetry off".into(),
                collocation_ms: "n/a".into(),
                mesh_refinement_ms: "n/a; fixed mesh".into(),
                solver_internal_ms: "n/a; telemetry off".into(),
                full_solve_ms: "n/a; use outer solve wall".into(),
                callback_parallel_dispatches: "n/a".into(),
                callback_sequential_dispatches: "n/a".into(),
                superlu_factor_calls: "n/a".into(),
                superlu_factor_ms: "n/a".into(),
                superlu_solve_ms: "n/a".into(),
            });
        }
        BvpSolverTelemetrySnapshot::Damped(snapshot) => {
            let c = &snapshot.counters;
            let t = &snapshot.timings;
            let callback = snapshot
                .atom_lambdify
                .as_ref()
                .or(snapshot.legacy_lambdify.as_ref());
            rows.push(NativeRow {
                scenario: scenario.name.into(),
                route: route.label(),
                execution: execution.into(),
                residual_or_rhs_calls: c.residual_requests.to_string(),
                jacobian_calls_or_builds: c.jacobian_requests.to_string(),
                linear_assemblies: "n/a".into(),
                factorizations: c.factorizations.to_string(),
                linear_solves: c.linear_solves.to_string(),
                collocation_evaluations: "n/a".into(),
                mesh_refinements: c.grid_refinements.to_string(),
                residual_or_rhs_ms: format!("{:.4}", ms(t.residual)),
                jacobian_ms: format!("{:.4}", ms(t.jacobian)),
                linear_assembly_ms: "n/a".into(),
                factorization_ms: format!("{:.4}", ms(t.factorization)),
                linear_solve_ms: format!("{:.4}", ms(t.linear_system)),
                collocation_ms: "n/a".into(),
                mesh_refinement_ms: format!("{:.4}", ms(t.grid_refinement)),
                solver_internal_ms: format!("{:.4}", ms(t.symbolic_operations)),
                full_solve_ms: "n/a; use outer solve wall".into(),
                callback_parallel_dispatches: callback.map_or_else(
                    || "n/a".into(),
                    |snapshot| snapshot.parallel_dispatches.to_string(),
                ),
                callback_sequential_dispatches: callback.map_or_else(
                    || "n/a".into(),
                    |snapshot| snapshot.sequential_dispatches.to_string(),
                ),
                superlu_factor_calls: "n/a".into(),
                superlu_factor_ms: "n/a".into(),
                superlu_solve_ms: "n/a".into(),
            });
        }
        BvpSolverTelemetrySnapshot::Sci(snapshot) => {
            rows.push(NativeRow {
                scenario: scenario.name.into(),
                route: route.label(),
                execution: execution.into(),
                residual_or_rhs_calls: snapshot.rhs_calls.to_string(),
                jacobian_calls_or_builds: snapshot.newton_jacobian_refreshes.to_string(),
                linear_assemblies: snapshot.linear_assemblies.to_string(),
                factorizations: snapshot.factorizations.to_string(),
                linear_solves: snapshot.linear_solves.to_string(),
                collocation_evaluations: snapshot.collocation_evaluations.to_string(),
                mesh_refinements: snapshot.mesh_refinements.to_string(),
                residual_or_rhs_ms: optional_ms(snapshot.rhs_callback_ms),
                jacobian_ms: optional_ms(snapshot.jacobian_callback_ms),
                linear_assembly_ms: optional_ms(snapshot.linear_assembly_ms),
                factorization_ms: optional_ms(snapshot.factorization_ms),
                linear_solve_ms: optional_ms(snapshot.linear_solve_ms),
                collocation_ms: optional_ms(snapshot.collocation_ms),
                mesh_refinement_ms: optional_ms(snapshot.mesh_refinement_ms),
                solver_internal_ms: optional_ms(snapshot.newton_ms),
                full_solve_ms: optional_ms(snapshot.full_solve_ms),
                callback_parallel_dispatches: snapshot.parallel_dispatches.to_string(),
                callback_sequential_dispatches: snapshot.sequential_dispatches.to_string(),
                superlu_factor_calls: snapshot.superlu.native_factorization_calls.to_string(),
                superlu_factor_ms: optional_ms(snapshot.superlu.native_factorization_ms),
                superlu_solve_ms: optional_ms(snapshot.superlu.native_solve_ms),
            });
            if snapshot.mesh_passes.is_empty() && snapshot.full_solve_ms.unwrap_or_default() > 0.0 {
                eprintln!(
                    "{} {} Sci telemetry had elapsed time but no mesh-pass detail",
                    scenario.name,
                    route.label()
                );
            }
        }
        BvpSolverTelemetrySnapshot::Frozen(snapshot) => {
            let c = &snapshot.counters;
            let t = &snapshot.timings;
            let callback = snapshot
                .atom_lambdify
                .as_ref()
                .or(snapshot.legacy_lambdify.as_ref());
            rows.push(NativeRow {
                scenario: scenario.name.into(),
                route: route.label(),
                execution: execution.into(),
                residual_or_rhs_calls: c.residual_requests.to_string(),
                jacobian_calls_or_builds: c.jacobian_requests.to_string(),
                linear_assemblies: "n/a".into(),
                factorizations: c.factorizations.to_string(),
                linear_solves: c.linear_solves.to_string(),
                collocation_evaluations: "n/a".into(),
                mesh_refinements: c.grid_refinements.to_string(),
                residual_or_rhs_ms: format!("{:.4}", ms(t.residual)),
                jacobian_ms: format!("{:.4}", ms(t.jacobian)),
                linear_assembly_ms: "n/a".into(),
                factorization_ms: format!("{:.4}", ms(t.factorization)),
                linear_solve_ms: format!("{:.4}", ms(t.linear_system)),
                collocation_ms: "n/a".into(),
                mesh_refinement_ms: format!("{:.4}", ms(t.grid_refinement)),
                solver_internal_ms: format!("{:.4}", ms(t.symbolic_operations)),
                full_solve_ms: format!("{:.4}", ms(t.total)),
                callback_parallel_dispatches: callback.map_or_else(
                    || "n/a".into(),
                    |snapshot| snapshot.parallel_dispatches.to_string(),
                ),
                callback_sequential_dispatches: callback.map_or_else(
                    || "n/a".into(),
                    |snapshot| snapshot.sequential_dispatches.to_string(),
                ),
                superlu_factor_calls: "n/a".into(),
                superlu_factor_ms: "n/a".into(),
                superlu_solve_ms: "n/a".into(),
            });
        }
    }
    rows
}

fn diagnostic_rows(
    scenario: Scenario,
    route: Route,
    execution: &str,
    samples: &[Sample],
) -> (MacroRow, Vec<StageRow>, NativeRow) {
    let mut equation = samples
        .iter()
        .map(|sample| ms(sample.equation_build))
        .collect::<Vec<_>>();
    let mut problem = samples
        .iter()
        .map(|sample| ms(sample.problem_validation))
        .collect::<Vec<_>>();
    let mut preparation = samples
        .iter()
        .map(|sample| ms(sample.solver_preparation))
        .collect::<Vec<_>>();
    let mut solve = samples
        .iter()
        .map(|sample| ms(sample.native_solve))
        .collect::<Vec<_>>();
    let mut total = samples
        .iter()
        .map(|sample| ms(sample.total))
        .collect::<Vec<_>>();
    let representative = &samples[samples.len() / 2];
    let preparation_ms = median(&mut preparation);
    let solve_ms = median(&mut solve);
    let total_ms = median(&mut total);
    let macro_row = MacroRow {
        scenario: scenario.name.into(),
        route: route.label(),
        execution: execution.into(),
        states: scenario.states,
        nodes: representative.mesh_nodes,
        max_error: format!(
            "{:.3e}",
            samples
                .iter()
                .map(|sample| sample.max_analytic_error)
                .fold(0.0, f64::max)
        ),
        preparation_ms: format!("{preparation_ms:.4}"),
        solve_ms: format!("{solve_ms:.4}"),
        total_ms: format!("{total_ms:.4}"),
    };
    let stages = [
        ("equation construction", median(&mut equation)),
        ("BvpProblem validation", median(&mut problem)),
        ("solver preparation", preparation_ms),
        ("native solve", solve_ms),
    ]
    .into_iter()
    .map(|(stage, time)| StageRow {
        scenario: scenario.name.into(),
        route: route.label(),
        execution: execution.into(),
        stage: stage.into(),
        median_ms: format!("{time:.4}"),
    })
    .collect();
    let native = collect_native_rows(scenario, route, execution, &representative.telemetry)
        .into_iter()
        .next()
        .expect("route telemetry is retained");
    (macro_row, stages, native)
}

fn callback_dispatches(telemetry: &BvpSolverTelemetrySnapshot) -> (u64, u64, f64) {
    match telemetry {
        BvpSolverTelemetrySnapshot::CorrectionNewton(_) => (0, 0, 0.0),
        BvpSolverTelemetrySnapshot::Damped(snapshot)
        | BvpSolverTelemetrySnapshot::Frozen(snapshot) => {
            let callback = snapshot
                .atom_lambdify
                .as_ref()
                .or(snapshot.legacy_lambdify.as_ref());
            callback.map_or((0, 0, 0.0), |callback| {
                (
                    callback.parallel_dispatches,
                    callback.sequential_dispatches,
                    (callback.residual_elapsed + callback.jacobian_elapsed).as_secs_f64() * 1000.0,
                )
            })
        }
        BvpSolverTelemetrySnapshot::Sci(snapshot) => (
            snapshot.parallel_dispatches,
            snapshot.sequential_dispatches,
            snapshot.rhs_callback_ms.unwrap_or_default()
                + snapshot.jacobian_callback_ms.unwrap_or_default(),
        ),
    }
}

fn large_callback_policy_rows() -> Vec<CallbackPolicyRow> {
    let scenario = *SCENARIOS.last().expect("large scenario exists");
    let routes = [
        Route::new(Family::Damped, Frontend::AtomView, Layout::Sparse),
        Route::new(Family::Sci, Frontend::AtomView, Layout::Sparse),
    ];
    let mut rows = Vec::new();
    for route in routes {
        let mut by_policy = Vec::new();
        for policy in [
            CallbackPolicy::Sequential,
            CallbackPolicy::Auto,
            CallbackPolicy::ForcedParallel,
        ] {
            let samples = (0..DIAGNOSTIC_REPETITIONS)
                .map(|_| run_once_with_policy(scenario, route, true, None, Some(policy)))
                .collect::<Vec<_>>();
            by_policy.push((policy, samples));
        }
        let reference = by_policy[0].1[by_policy[0].1.len() / 2].values.clone();
        for (policy, samples) in by_policy {
            let representative = &samples[samples.len() / 2];
            let max_profile_delta = samples
                .iter()
                .flat_map(|sample| sample.values.iter().zip(&reference))
                .map(|(actual, expected)| (actual - expected).abs())
                .fold(0.0, f64::max);
            assert!(
                max_profile_delta <= 2e-5,
                "{} callback policy {} changed the profile by {max_profile_delta:e}",
                route.label(),
                policy.label()
            );
            let (parallel, sequential, callback_ms) =
                callback_dispatches(&representative.telemetry);
            if matches!(policy, CallbackPolicy::ForcedParallel) {
                assert!(
                    parallel > 0,
                    "forced parallel policy did not dispatch work for {}",
                    route.label()
                );
            }
            let mut solve_times = samples
                .iter()
                .map(|sample| ms(sample.native_solve))
                .collect::<Vec<_>>();
            rows.push(CallbackPolicyRow {
                scenario: scenario.name.into(),
                solver: route.family.label().into(),
                policy: policy.label().into(),
                max_profile_delta: format!("{max_profile_delta:.3e}"),
                parallel_dispatches: parallel,
                sequential_dispatches: sequential,
                callback_ms: format!("{callback_ms:.4}"),
                median_solve_ms: format!("{:.4}", median(&mut solve_times)),
            });
        }
    }
    rows
}

fn report_detailed_axes() {
    let mut macros = Vec::new();
    let mut stages = Vec::new();
    let mut native = Vec::new();
    let mut aot_comparisons = Vec::new();

    for scenario in SCENARIOS {
        for route in routes() {
            let lambda_samples = (0..DIAGNOSTIC_REPETITIONS)
                .map(|_| run_once(scenario, route, true, None))
                .collect::<Vec<_>>();
            let (macro_row, stage_rows, native_row) =
                diagnostic_rows(scenario, route, "Lambdify", &lambda_samples);
            let lambda_total = ms(lambda_samples[lambda_samples.len() / 2].total);
            macros.push(macro_row);
            stages.extend(stage_rows);
            native.push(native_row);

            if is_aot_enabled() && supports_aot(route) {
                let output = tempfile::tempdir().expect("route AOT artifact directory");
                let aot_samples = (0..DIAGNOSTIC_REPETITIONS)
                    .map(|_| run_once(scenario, route, true, Some(output.path())))
                    .collect::<Vec<_>>();
                let lambda_values = &lambda_samples[lambda_samples.len() / 2].values;
                for sample in &aot_samples {
                    assert_eq!(sample.values.len(), lambda_values.len());
                    let max_delta = sample
                        .values
                        .iter()
                        .zip(lambda_values)
                        .map(|(aot, lambda)| (aot - lambda).abs())
                        .fold(0.0, f64::max);
                    assert!(
                        max_delta <= 2e-5,
                        "AOT profile differs from Lambdify for {} / {}: max delta {max_delta:e}",
                        scenario.name,
                        route.label()
                    );
                }
                let (macro_row, stage_rows, native_row) =
                    diagnostic_rows(scenario, route, "AOT-TCC", &aot_samples);
                let aot_total = ms(aot_samples[aot_samples.len() / 2].total);
                aot_comparisons.push(ComparisonRow {
                    scenario: scenario.name.into(),
                    route: route.label(),
                    lambdify_total_ms: format!("{lambda_total:.4}"),
                    aot_total_ms: format!("{aot_total:.4}"),
                    aot_delta_pct: format!("{:+.1}", 100.0 * (aot_total / lambda_total - 1.0)),
                });
                macros.push(macro_row);
                stages.extend(stage_rows);
                native.push(native_row);
            }
        }
    }

    println!(
        "\nPERF-BVP-01 | analytic profiles y_i'=-k_i*y_i, y_i(0)=1, y_i in [0,1] | cases: scalar/24, 8-state/64, 24-state/128"
    );
    println!(
        "\nDetailed route summary (three Detailed repetitions; these are attribution runs, not Criterion wall samples)\n{}",
        Table::new(macros).with(Style::markdown())
    );
    println!(
        "\nStage attribution: equation parsing, API problem validation, native solver preparation, native solve\n{}",
        Table::new(stages).with(Style::markdown())
    );
    println!(
        "\nNative BVP telemetry (times are route-native; Damped and Sci timer scopes have different semantics)\n{}",
        Table::new(native).with(Style::markdown())
    );
    println!(
        "\nKiThe-side preparation before the RST handoff (40 species / 39 reactions; equivalent equations asserted across policies)\n{}",
        Table::new(reactor_setup_policy_rows()).with(Style::markdown())
    );
    println!(
        "\nLargest-case callback policy comparison (AtomView/Sparse; timing and dispatches verified independently for Damped and BVP_sci)\n{}",
        Table::new(large_callback_policy_rows()).with(Style::markdown())
    );
    if !aot_comparisons.is_empty() {
        println!(
            "\nAOT-TCC versus Lambdify, matched solver/frontend/layout/scenario; delta is instrumented end-to-end total\n{}",
            Table::new(aot_comparisons).with(Style::markdown())
        );
    } else {
        println!("\nAOT axis skipped: set KITHE_RUN_BVP_AOT_BENCH=1 and install tcc.");
    }
}

fn reactor_bvp_baseline(c: &mut Criterion) {
    report_detailed_axes();
    let mut group = c.benchmark_group("reactor_bvp_analytic_profiles");
    group.sample_size(25);
    for scenario in SCENARIOS {
        for route in routes() {
            let id = format!("{}/{}", scenario.name, route.label());
            group.bench_function(BenchmarkId::new("lambdify", id), |b| {
                b.iter(|| {
                    let sample = run_once(scenario, route, false, None);
                    black_box(sample.values.last().copied().unwrap_or_default())
                });
            });
        }
    }
    let scenario = *SCENARIOS.last().expect("large scenario exists");
    for route in [
        Route::new(Family::Damped, Frontend::AtomView, Layout::Sparse),
        Route::new(Family::Sci, Frontend::AtomView, Layout::Sparse),
    ] {
        for policy in [
            CallbackPolicy::Sequential,
            CallbackPolicy::Auto,
            CallbackPolicy::ForcedParallel,
        ] {
            let id = format!("{}/{}", route.label(), policy.label());
            group.bench_function(BenchmarkId::new("callback-policy", id), |b| {
                b.iter(|| {
                    let sample = run_once_with_policy(scenario, route, false, None, Some(policy));
                    black_box(sample.values.last().copied().unwrap_or_default())
                });
            });
        }
    }
    group.finish();

    let mut setup_group = c.benchmark_group("reactor_bvp_kithe_preparation");
    setup_group.sample_size(25);
    for policy in [
        SymbolicRhsAssemblyPolicy::Sequential,
        SymbolicRhsAssemblyPolicy::Auto,
        SymbolicRhsAssemblyPolicy::Parallel,
    ] {
        setup_group.bench_function(
            BenchmarkId::new("40species-39reactions", format!("{policy:?}")),
            |b| {
                b.iter(|| {
                    let mut reactor = large_reactor_setup_fixture(policy);
                    reactor
                        .setup_bvp()
                        .expect("large KiThe-side BVP setup should complete");
                    black_box(reactor.solver.eq_system.len())
                });
            },
        );
    }
    setup_group.finish();

    let mut physical_group = c.benchmark_group("reactor_bvp_large_mechanism_physical");
    physical_group.sample_size(10);
    physical_group.measurement_time(Duration::from_secs(12));
    for (label, method) in [
        ("DoublePoints", GridRefinementMethod::DoublePoints),
        ("Pearson-1e-2", GridRefinementMethod::Pearson(1e-2, 1.4)),
        (
            "GrcarSmooke-1e-2-5e-2",
            GridRefinementMethod::GrcarSmooke(1e-2, 5e-2, 1.4),
        ),
        (
            "TwoPoint-2e-2-C1.4",
            GridRefinementMethod::TwoPoint(2e-2, 2e-2, 1.4),
        ),
    ] {
        physical_group.bench_function(BenchmarkId::new("Damped", label), |b| {
            b.iter_batched(
                large_mechanism_fixture::build_reactor,
                |mut reactor| {
                    let node_count = solve_large_mechanism_damped(&mut reactor, method.clone());
                    black_box(node_count)
                },
                BatchSize::SmallInput,
            );
        });
    }
    physical_group.bench_function(BenchmarkId::new("CorrectionNewton", "fixed-65"), |b| {
        b.iter_batched(
            large_mechanism_fixture::build_reactor,
            |mut reactor| {
                let node_count = solve_large_mechanism_correction_newton(&mut reactor, 65);
                black_box(node_count)
            },
            BatchSize::SmallInput,
        );
    });
    physical_group.finish();
}

fn large_mechanism_problem(
    reactor: &mut KiThe::SimpleReactorsBVP::SimpleReactorBVP::SimpleReactorTask,
    interval_count: usize,
) -> BvpProblem {
    reactor
        .setup_bvp()
        .expect("large-mechanism KiThe BVP preparation");
    reactor.solver.arg_name = "z".to_string();
    reactor.solver.x_range = (0.0, 1.0);
    let mesh = (0..=interval_count)
        .map(|index| index as f64 / interval_count as f64)
        .collect::<Vec<_>>();
    let conditions = reactor
        .solver
        .unknowns
        .iter()
        .enumerate()
        .map(|(state, name)| {
            let (side, value) = reactor.solver.BorderConditions[name];
            BvpBoundaryCondition {
                state,
                side: if side == 0 {
                    BvpBoundarySide::Left
                } else {
                    BvpBoundarySide::Right
                },
                value,
            }
        })
        .collect();
    BvpProblem::new(
        reactor.solver.eq_system.clone(),
        reactor.solver.unknowns.clone(),
        "z",
        mesh,
        conditions,
    )
    .expect("large-mechanism BVP problem validates")
}

fn solve_large_mechanism_damped(
    reactor: &mut KiThe::SimpleReactorsBVP::SimpleReactorBVP::SimpleReactorTask,
    method: GridRefinementMethod,
) -> usize {
    let interval_count = 33;
    let problem = large_mechanism_problem(reactor, interval_count);
    let initial_guess = large_mechanism_fixture::inlet_initial_guess(reactor, interval_count);
    let mut options = ReactorBvpSolverConfig::sparse_lambdify()
        .with_symbolic_backend(ReactorBvpSymbolicBackend::ExprLegacy)
        .to_rusted_options()
        .expect("Damped symbolic options");
    options = options
        .with_scheme_name("trapezoid".to_string())
        .with_strategy_params(Some(SolverParams {
            max_jac: Some(1),
            max_damp_iter: Some(12),
            adaptive: Some(AdaptiveGridConfig {
                version: 1,
                max_refinements: 2,
                grid_method: method,
            }),
            ..SolverParams::default()
        }))
        .with_abs_tolerance(1e-5)
        .with_max_iterations(150)
        .with_bvp_telemetry_mode(BvpTelemetryMode::Off);
    let relative_tolerances = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| (name.clone(), 1e-5))
        .collect();
    let bounds = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| {
            let bound = if name == "Teta" {
                (-1.0, 3.0)
            } else if name == "q" || name.starts_with('J') {
                (-1e4, 1e4)
            } else {
                (-0.1, 1.1)
            };
            (name.clone(), bound)
        })
        .collect();
    options = options
        .with_rel_tolerance(relative_tolerances)
        .with_bounds(bounds);
    let config = BvpSolverConfig::Damped {
        options,
        initial_guess,
    };
    let mut solver = BvpSolver::new(problem, config).expect("Damped preparation");
    solver.solve().expect("physical Damped refinement solve");
    let view = solver.solution_view().expect("Damped solution view");
    assert!(view.mesh().len() >= interval_count + 1);
    for node in 0..view.mesh().len() {
        for state in 0..reactor.solver.unknowns.len() {
            assert!(
                view.value(node, state).is_some_and(f64::is_finite),
                "Damped result contains a non-finite state"
            );
        }
    }
    view.mesh().len()
}

fn solve_large_mechanism_correction_newton(
    reactor: &mut KiThe::SimpleReactorsBVP::SimpleReactorBVP::SimpleReactorTask,
    interval_count: usize,
) -> usize {
    let problem = large_mechanism_problem(reactor, interval_count);
    let initial_guess = large_mechanism_fixture::inlet_initial_guess(reactor, interval_count);
    let mut options = CorrectionNewtonOptions {
        symbolic_frontend: CorrectionNewtonSymbolicFrontend::ExprLegacy,
        matrix_layout: CorrectionNewtonMatrixLayout::Sparse,
        absolute_tolerance: 1e-5,
        default_relative_tolerance: 1e-5,
        max_iterations: 100,
        telemetry_mode: CorrectionNewtonTelemetryMode::Off,
        ..CorrectionNewtonOptions::default()
    };
    options.relative_tolerances = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| (name.clone(), 1e-5))
        .collect();
    options.bounds = reactor
        .solver
        .unknowns
        .iter()
        .map(|name| {
            let bound = if name == "Teta" {
                (-1.0, 3.0)
            } else if name == "q" || name.starts_with('J') {
                (-1e4, 1e4)
            } else {
                (0.0, 1.1)
            };
            (name.clone(), bound)
        })
        .collect();
    let config = BvpSolverConfig::CorrectionNewton {
        options,
        initial_guess,
    };
    let mut solver = BvpSolver::new(problem, config).expect("CorrectionNewton preparation");
    solver
        .solve()
        .expect("physical CorrectionNewton fixed-mesh solve");
    let view = solver
        .solution_view()
        .expect("CorrectionNewton solution view");
    assert!(view.mesh().len() >= interval_count + 1);
    for node in 0..view.mesh().len() {
        for state in 0..reactor.solver.unknowns.len() {
            assert!(
                view.value(node, state).is_some_and(f64::is_finite),
                "CorrectionNewton result contains a non-finite state"
            );
        }
    }
    view.mesh().len()
}

criterion_group!(benches, reactor_bvp_baseline);
criterion_main!(benches);
