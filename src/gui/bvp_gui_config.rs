use crate::ReactorsBVP::solver_backend::{
    ReactorBvpAotBuildPolicy, ReactorBvpAotBuildProfile, ReactorBvpAotCompiler,
    ReactorBvpAotConfig, ReactorBvpAotExecutionPolicy, ReactorBvpExecutionBackend,
    ReactorBvpMatrixBackend, ReactorBvpSolverConfig, ReactorBvpSymbolicBackend,
};
use crate::ReactorsBVP::task_value_conversion::try_value_as_usize;
use RustedSciThe::command_interpreter::task_parser::{DocumentMap, Value};
use log::info;
use std::collections::HashMap;

/// Migration notes produced while folding a raw document into the typed GUI model.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct BvpGuiMigrationReport {
    pub warnings: Vec<String>,
}

impl BvpGuiMigrationReport {
    fn push_warning(&mut self, message: impl Into<String>) {
        self.warnings.push(message.into());
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum BvpScheme {
    #[default]
    Forward,
    Trapezoid,
}

impl BvpScheme {
    fn parse(raw: &str) -> Option<Self> {
        match normalize_key(raw).as_str() {
            "forward" => Some(Self::Forward),
            "trapezoid" => Some(Self::Trapezoid),
            _ => None,
        }
    }

    fn as_document_value(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Trapezoid => "trapezoid",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum BvpStrategy {
    #[default]
    Damped,
    Frozen,
    Naive,
}

impl BvpStrategy {
    fn parse(raw: &str) -> Option<Self> {
        match normalize_key(raw).as_str() {
            "damped" => Some(Self::Damped),
            "frozen" => Some(Self::Frozen),
            "naive" => Some(Self::Naive),
            _ => None,
        }
    }

    fn as_document_value(self) -> &'static str {
        match self {
            Self::Damped => "Damped",
            Self::Frozen => "Frozen",
            Self::Naive => "Naive",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum BvpBandedLinearSolver {
    #[default]
    Auto,
    Faithful,
    BlockTridiagonal,
    BlockTridiagonalConsistent,
    FaerSparse,
    GeneralPartialPivot,
}

impl BvpBandedLinearSolver {
    fn parse(raw: &str) -> Option<Self> {
        match normalize_key(raw).as_str() {
            "auto" => Some(Self::Auto),
            "faithful" => Some(Self::Faithful),
            "block_tridiagonal" => Some(Self::BlockTridiagonal),
            "block_tridiagonal_consistent" => Some(Self::BlockTridiagonalConsistent),
            "faer_sparse" => Some(Self::FaerSparse),
            "general_partial_pivot" => Some(Self::GeneralPartialPivot),
            _ => None,
        }
    }

    fn as_document_value(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Faithful => "faithful",
            Self::BlockTridiagonal => "block_tridiagonal",
            Self::BlockTridiagonalConsistent => "block_tridiagonal_consistent",
            Self::FaerSparse => "faer_sparse",
            Self::GeneralPartialPivot => "general_partial_pivot",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum BvpGridRefinementMethod {
    #[default]
    Pearson,
    GrcarSmooke,
    Twopnt,
    Easy,
    DoublePoints,
}

impl BvpGridRefinementMethod {
    fn parse(raw: &str) -> Option<Self> {
        match normalize_key(raw).as_str() {
            "pearson" => Some(Self::Pearson),
            "grcar_smooke" | "grcarsmooke" => Some(Self::GrcarSmooke),
            "twopnt" => Some(Self::Twopnt),
            "easy" => Some(Self::Easy),
            "double_points" | "doublepoints" => Some(Self::DoublePoints),
            _ => None,
        }
    }

    fn as_document_key(self) -> &'static str {
        match self {
            Self::Pearson => "pearson",
            Self::GrcarSmooke => "grcarsmooke",
            Self::Twopnt => "twopnt",
            Self::Easy => "easy",
            Self::DoublePoints => "doubleoints",
        }
    }

    fn default_params(self) -> Vec<f64> {
        match self {
            Self::Pearson => vec![0.05, 2.5],
            Self::GrcarSmooke => vec![0.05, 0.05, 1.25],
            Self::Twopnt => vec![0.05, 0.05, 1.25],
            Self::Easy => vec![0.05],
            Self::DoublePoints => Vec::new(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpGridRefinementConfig {
    pub method: BvpGridRefinementMethod,
    pub params: Vec<f64>,
}

impl Default for BvpGridRefinementConfig {
    fn default() -> Self {
        let method = BvpGridRefinementMethod::default();
        Self {
            params: method.default_params(),
            method,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpAdaptiveGridConfig {
    pub enabled: bool,
    pub version: usize,
    pub max_refinements: usize,
    pub grid_refinement: BvpGridRefinementConfig,
}

impl Default for BvpAdaptiveGridConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            version: 1,
            max_refinements: 3,
            grid_refinement: BvpGridRefinementConfig::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum BvpInitialGuessMode {
    #[default]
    Universal,
    PerVariable,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpInitialGuessConfig {
    pub mode: BvpInitialGuessMode,
    pub universal: f64,
    pub c: f64,
    pub j: f64,
    pub teta: f64,
    pub q: f64,
}

impl Default for BvpInitialGuessConfig {
    fn default() -> Self {
        Self {
            mode: BvpInitialGuessMode::Universal,
            universal: 1e-2,
            c: 1e-2,
            j: 1e-2,
            teta: 1e-2,
            q: 1e-2,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpRelativeToleranceConfig {
    pub c: f64,
    pub j: f64,
    pub teta: f64,
    pub q: f64,
}

impl Default for BvpRelativeToleranceConfig {
    fn default() -> Self {
        Self {
            c: 1e-7,
            j: 1e-7,
            teta: 1e-7,
            q: 1e-7,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpStrategyParamsConfig {
    pub max_jac: Option<usize>,
    pub max_damp_iter: Option<usize>,
    pub damp_factor: Option<f64>,
}

impl Default for BvpStrategyParamsConfig {
    fn default() -> Self {
        Self {
            max_jac: Some(3),
            max_damp_iter: Some(10),
            damp_factor: Some(0.5),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpAdvancedSolverConfig {
    pub rel_tolerance: BvpRelativeToleranceConfig,
    pub strategy_params: BvpStrategyParamsConfig,
}

impl Default for BvpAdvancedSolverConfig {
    fn default() -> Self {
        Self {
            rel_tolerance: BvpRelativeToleranceConfig::default(),
            strategy_params: BvpStrategyParamsConfig::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpCoreSolverConfig {
    pub scheme: BvpScheme,
    pub strategy: BvpStrategy,
    pub linear_sys_method: Option<String>,
    pub abs_tolerance: f64,
    pub max_iterations: usize,
    pub loglevel: Option<String>,
    pub dont_save_log: bool,
    pub banded_linear_solver: BvpBandedLinearSolver,
    pub refinement_steps: usize,
}

impl Default for BvpCoreSolverConfig {
    fn default() -> Self {
        Self {
            scheme: BvpScheme::Forward,
            strategy: BvpStrategy::Damped,
            linear_sys_method: None,
            abs_tolerance: 1e-7,
            max_iterations: 100,
            loglevel: Some("info".to_string()),
            dont_save_log: true,
            banded_linear_solver: BvpBandedLinearSolver::Auto,
            refinement_steps: 5,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpPostprocessingConfig {
    pub plot: bool,
    pub gnuplot: bool,
    pub save_to_csv: bool,
    pub filename: String,
    pub save: bool,
    pub return_to_dimension: bool,
    pub no_plots_in_terminal: bool,
    pub gui_plot: bool,
}

impl Default for BvpPostprocessingConfig {
    fn default() -> Self {
        Self {
            plot: true,
            gnuplot: false,
            save_to_csv: false,
            filename: String::new(),
            save: false,
            return_to_dimension: false,
            no_plots_in_terminal: true,
            gui_plot: true,
        }
    }
}

/// Typed snapshot of the BVP GUI structured state.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BvpGuiConfig {
    pub core: BvpCoreSolverConfig,
    pub backend: ReactorBvpSolverConfig,
    pub advanced_solver: BvpAdvancedSolverConfig,
    pub initial_guess: BvpInitialGuessConfig,
    pub adaptive_grid: BvpAdaptiveGridConfig,
    pub postprocessing: BvpPostprocessingConfig,
}

impl Default for BvpGuiConfig {
    fn default() -> Self {
        Self {
            core: BvpCoreSolverConfig::default(),
            backend: ReactorBvpSolverConfig::default_lambdify(),
            advanced_solver: BvpAdvancedSolverConfig::default(),
            initial_guess: BvpInitialGuessConfig::default(),
            adaptive_grid: BvpAdaptiveGridConfig::default(),
            postprocessing: BvpPostprocessingConfig::default(),
        }
    }
}

impl BvpGuiConfig {
    /// Fold a raw task document into the structured GUI model.
    pub(crate) fn from_document(document: &DocumentMap) -> (Self, BvpGuiMigrationReport) {
        let mut report = BvpGuiMigrationReport::default();
        let solver_section = document.get("solver_settings");
        let postprocessing = document.get("postprocessing");
        let adaptive_strategy = document.get("adaptive_strategy");
        let grid_refinement = document.get("grid_refinement");

        let mut config = Self::default();
        if let Some(section) = solver_section {
            config.core = parse_core_solver_config(section, &mut report);
            config.backend = parse_backend_config(section, &mut report);
        }
        config.advanced_solver = parse_advanced_solver_config(
            document.get("rel_tolerance"),
            document.get("strategy_params"),
            &mut report,
        );
        config.initial_guess =
            parse_initial_guess_config(document.get("initial_guess"), &mut report);
        config.postprocessing =
            parse_postprocessing_config(postprocessing, solver_section, &mut report);
        config.adaptive_grid =
            parse_adaptive_grid_config(document, adaptive_strategy, grid_refinement, &mut report);

        (config, report)
    }

    /// Write the typed GUI model back into the raw task document.
    pub(crate) fn apply_to_document(&self, document: &mut DocumentMap) {
        {
            let solver_section = ensure_section_mut(document, "solver_settings");
            self.core.apply_to_section(solver_section);
            self.backend.apply_to_section(solver_section);
            for field_name in BVP_POSTPROCESSING_FIELD_NAMES {
                solver_section.remove(*field_name);
            }
        }

        self.advanced_solver.apply_to_document(document);
        self.initial_guess.apply_to_document(document);

        let postprocessing_section = ensure_section_mut(document, "postprocessing");
        self.postprocessing.apply_to_section(postprocessing_section);
        self.adaptive_grid.apply_to_document(document);
    }
}

const BVP_POSTPROCESSING_FIELD_NAMES: &[&str] = &[
    "plot",
    "gnuplot",
    "save_to_csv",
    "filename",
    "save",
    "return_to_dimension",
    "no_plots_in_terminal",
    "gui_plot",
];

impl ReactorBvpSolverConfig {
    fn apply_to_section(&self, section: &mut HashMap<String, Option<Vec<Value>>>) {
        set_string_field(section, "generated_backend", &generated_backend_name(self));
        set_string_field(
            section,
            "matrix_backend",
            match self.matrix_backend {
                ReactorBvpMatrixBackend::Banded => "Banded",
                ReactorBvpMatrixBackend::Sparse => "Sparse",
            },
        );
        set_string_field(
            section,
            "method",
            match self.matrix_backend {
                ReactorBvpMatrixBackend::Banded => "Banded",
                ReactorBvpMatrixBackend::Sparse => "Sparse",
            },
        );
        set_string_field(
            section,
            "symbolic_backend",
            match self.symbolic_backend {
                ReactorBvpSymbolicBackend::AtomView => "AtomView",
                ReactorBvpSymbolicBackend::ExprLegacy => "ExprLegacy",
            },
        );

        match self.execution_backend {
            ReactorBvpExecutionBackend::Lambdify => {
                set_string_field(section, "backend_policy", "lambdify_only");
                section.remove("aot_c_compiler");
                section.remove("aot_build_policy");
                section.remove("aot_build_profile");
                section.remove("aot_compile_preset");
                section.remove("aot_execution_policy");
                section.remove("aot_residual_target_chunks");
                section.remove("aot_jacobian_target_chunks");
            }
            ReactorBvpExecutionBackend::Aot(aot) => {
                set_string_field(section, "backend_policy", "aot_only");
                set_string_field(
                    section,
                    "aot_c_compiler",
                    match aot.compiler {
                        ReactorBvpAotCompiler::CTcc => "tcc",
                        ReactorBvpAotCompiler::CGcc => "gcc",
                        ReactorBvpAotCompiler::Zig => "zig",
                    },
                );
                set_string_field(
                    section,
                    "aot_build_policy",
                    match aot.build_policy {
                        ReactorBvpAotBuildPolicy::UseIfAvailable => "use_if_available",
                        ReactorBvpAotBuildPolicy::BuildIfMissing => "build_if_missing",
                        ReactorBvpAotBuildPolicy::RequirePrebuilt => "require_prebuilt",
                        ReactorBvpAotBuildPolicy::RebuildAlways => "rebuild_always",
                    },
                );
                set_string_field(
                    section,
                    "aot_build_profile",
                    match aot.build_profile {
                        ReactorBvpAotBuildProfile::Release => "release",
                        ReactorBvpAotBuildProfile::Debug => "debug",
                    },
                );
                set_string_field(
                    section,
                    "aot_compile_preset",
                    match aot.compile_preset {
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::Production =>
                            "production",
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::FastBuild =>
                            "fast_build",
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::DevFastest =>
                            "dev_fastest",
                    },
                );
                set_string_field(
                    section,
                    "aot_execution_policy",
                    match aot.execution_policy {
                        ReactorBvpAotExecutionPolicy::Auto => "auto",
                        ReactorBvpAotExecutionPolicy::Sequential => "sequential",
                        ReactorBvpAotExecutionPolicy::Parallel => "parallel",
                    },
                );
                if let Some(chunks) = aot.residual_target_chunks {
                    set_usize_field(section, "aot_residual_target_chunks", chunks);
                } else {
                    section.remove("aot_residual_target_chunks");
                }
                if let Some(chunks) = aot.jacobian_target_chunks {
                    set_usize_field(section, "aot_jacobian_target_chunks", chunks);
                } else {
                    section.remove("aot_jacobian_target_chunks");
                }
            }
        }
    }
}

impl BvpCoreSolverConfig {
    fn apply_to_section(&self, section: &mut HashMap<String, Option<Vec<Value>>>) {
        set_string_field(section, "scheme", self.scheme.as_document_value());
        set_string_field(section, "strategy", self.strategy.as_document_value());
        set_optional_string_field(
            section,
            "linear_sys_method",
            self.linear_sys_method.as_deref(),
        );
        set_float_field(section, "abs_tolerance", self.abs_tolerance);
        set_usize_field(section, "max_iterations", self.max_iterations);
        set_optional_string_field(section, "loglevel", self.loglevel.as_deref());
        set_bool_field(section, "dont_save_log", self.dont_save_log);
        set_string_field(
            section,
            "banded_linear_solver",
            self.banded_linear_solver.as_document_value(),
        );
        set_usize_field(section, "refinement_steps", self.refinement_steps);
    }
}

impl BvpPostprocessingConfig {
    fn apply_to_section(&self, section: &mut HashMap<String, Option<Vec<Value>>>) {
        set_bool_field(section, "plot", self.plot);
        set_bool_field(section, "gnuplot", self.gnuplot);
        set_bool_field(section, "save_to_csv", self.save_to_csv);
        set_string_field(section, "filename", &self.filename);
        set_bool_field(section, "save", self.save);
        set_bool_field(section, "return_to_dimension", self.return_to_dimension);
        set_bool_field(section, "no_plots_in_terminal", self.no_plots_in_terminal);
        set_bool_field(section, "gui_plot", self.gui_plot);
    }
}

impl BvpRelativeToleranceConfig {
    fn apply_to_section(&self, section: &mut HashMap<String, Option<Vec<Value>>>) {
        set_float_field(section, "C", self.c);
        set_float_field(section, "J", self.j);
        set_float_field(section, "Teta", self.teta);
        set_float_field(section, "q", self.q);
    }
}

impl BvpStrategyParamsConfig {
    fn apply_to_section(&self, section: &mut HashMap<String, Option<Vec<Value>>>) {
        set_optional_usize_field(section, "max_jac", self.max_jac);
        set_optional_usize_field(section, "max_damp_iter", self.max_damp_iter);
        set_optional_float_field(section, "damp_factor", self.damp_factor);
    }
}

impl BvpAdvancedSolverConfig {
    fn apply_to_document(&self, document: &mut DocumentMap) {
        let rel_tolerance = ensure_section_mut(document, "rel_tolerance");
        self.rel_tolerance.apply_to_section(rel_tolerance);

        let strategy_params = ensure_section_mut(document, "strategy_params");
        self.strategy_params.apply_to_section(strategy_params);
    }
}

impl BvpInitialGuessConfig {
    fn apply_to_document(&self, document: &mut DocumentMap) {
        let section = ensure_section_mut(document, "initial_guess");
        match self.mode {
            BvpInitialGuessMode::Universal => {
                set_float_field(section, "universal", self.universal);
                section.remove("C");
                section.remove("J");
                section.remove("Teta");
                section.remove("q");
            }
            BvpInitialGuessMode::PerVariable => {
                section.remove("universal");
                set_float_field(section, "C", self.c);
                set_float_field(section, "J", self.j);
                set_float_field(section, "Teta", self.teta);
                set_float_field(section, "q", self.q);
            }
        }
    }
}

impl BvpAdaptiveGridConfig {
    fn apply_to_document(&self, document: &mut DocumentMap) {
        if let Some(strategy_params) = document.get_mut("strategy_params") {
            strategy_params.remove("adaptive");
        }

        if !self.enabled {
            document.remove("adaptive_strategy");
            document.remove("grid_refinement");
            return;
        }

        let adaptive_strategy = ensure_section_mut(document, "adaptive_strategy");
        adaptive_strategy.insert(
            "version".to_string(),
            Some(vec![Value::Usize(self.version)]),
        );
        set_usize_field(adaptive_strategy, "max_refinements", self.max_refinements);

        let grid_refinement = ensure_section_mut(document, "grid_refinement");
        grid_refinement.clear();
        grid_refinement.insert(
            self.grid_refinement.method.as_document_key().to_string(),
            Some(vec![Value::Vector(self.grid_refinement.params.clone())]),
        );
    }
}

fn parse_core_solver_config(
    section: &HashMap<String, Option<Vec<Value>>>,
    report: &mut BvpGuiMigrationReport,
) -> BvpCoreSolverConfig {
    let mut core = BvpCoreSolverConfig::default();
    if let Some(scheme) = string_field(section, "scheme") {
        if let Some(parsed) = BvpScheme::parse(&scheme) {
            core.scheme = parsed;
        } else {
            report.push_warning(format!(
                "Unknown solver scheme `{scheme}`. Falling back to `forward`."
            ));
        }
    }
    if let Some(strategy) = string_field(section, "strategy") {
        if let Some(parsed) = BvpStrategy::parse(&strategy) {
            core.strategy = parsed;
        } else {
            report.push_warning(format!(
                "Unknown solver strategy `{strategy}`. Falling back to `Damped`."
            ));
        }
    }
    core.linear_sys_method = optional_string_field(section, "linear_sys_method");
    if let Some(value) = float_field(section, "abs_tolerance") {
        if value.is_finite() {
            core.abs_tolerance = value;
        } else {
            report.push_warning("Non-finite abs_tolerance encountered. Falling back to 1e-7.");
        }
    }
    if let Some(value) = usize_field(section, "max_iterations") {
        core.max_iterations = value;
    }
    core.loglevel = optional_string_field(section, "loglevel");
    if let Some(value) = bool_field(section, "dont_save_log") {
        core.dont_save_log = value;
    }
    if let Some(banded_solver) = string_field(section, "banded_linear_solver") {
        if let Some(parsed) = BvpBandedLinearSolver::parse(&banded_solver) {
            core.banded_linear_solver = parsed;
        } else {
            report.push_warning(format!(
                "Unknown banded linear solver `{banded_solver}`. Falling back to `auto`."
            ));
        }
    }
    if let Some(value) = usize_field(section, "refinement_steps") {
        core.refinement_steps = value;
    }

    core
}

fn parse_backend_config(
    section: &HashMap<String, Option<Vec<Value>>>,
    report: &mut BvpGuiMigrationReport,
) -> ReactorBvpSolverConfig {
    let backend_name = string_field(section, "generated_backend")
        .or_else(|| match string_field(section, "backend_policy") {
            Some(policy) if policy.trim().eq_ignore_ascii_case("aot_only") => {
                Some("banded_aot_tcc".to_string())
            }
            _ => None,
        })
        .unwrap_or_else(|| "banded_lambdify".to_string());

    let mut backend = match ReactorBvpSolverConfig::from_generated_backend_name(&backend_name) {
        Ok(config) => config,
        Err(error) => {
            report.push_warning(error);
            ReactorBvpSolverConfig::default_lambdify()
        }
    };

    if let Some(matrix_backend) =
        string_field(section, "matrix_backend").or_else(|| string_field(section, "method"))
    {
        match ReactorBvpMatrixBackend::parse(&matrix_backend) {
            Ok(parsed) => backend = backend.with_matrix_backend(parsed),
            Err(error) => report.push_warning(error),
        }
    }

    if let Some(symbolic_backend) = string_field(section, "symbolic_backend") {
        match ReactorBvpSymbolicBackend::parse(&symbolic_backend) {
            Ok(parsed) => backend = backend.with_symbolic_backend(parsed),
            Err(error) => report.push_warning(error),
        }
    }

    if matches!(
        backend.execution_backend,
        ReactorBvpExecutionBackend::Aot(_)
    ) {
        let mut aot_config = match backend.execution_backend {
            ReactorBvpExecutionBackend::Aot(config) => config,
            ReactorBvpExecutionBackend::Lambdify => ReactorBvpAotConfig::default(),
        };

        if let Some(compiler) = string_field(section, "aot_c_compiler") {
            match ReactorBvpAotCompiler::parse(&compiler) {
                Ok(parsed) => aot_config.compiler = parsed,
                Err(error) => report.push_warning(error),
            }
        }

        if let Some(build_policy) = string_field(section, "aot_build_policy") {
            match ReactorBvpAotBuildPolicy::parse(&build_policy) {
                Ok(parsed) => aot_config.build_policy = parsed,
                Err(error) => report.push_warning(error),
            }
        }
        if let Some(build_profile) = string_field(section, "aot_build_profile") {
            match ReactorBvpAotBuildProfile::parse(&build_profile) {
                Ok(parsed) => aot_config.build_profile = parsed,
                Err(error) => report.push_warning(error),
            }
        }
        if let Some(compile_preset) = string_field(section, "aot_compile_preset") {
            match normalize_key(&compile_preset).as_str() {
                "production" | "prod" => {
                    aot_config.compile_preset =
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::Production
                }
                "fast_build" | "fastbuild" => {
                    aot_config.compile_preset =
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::FastBuild
                }
                "dev_fastest" | "devfastest" | "fastest" => {
                    aot_config.compile_preset =
                        crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::DevFastest
                }
                other => report.push_warning(format!(
                    "Unknown AOT compile preset `{other}`. Falling back to `dev_fastest`."
                )),
            }
        }
        if let Some(execution_policy) = string_field(section, "aot_execution_policy") {
            match ReactorBvpAotExecutionPolicy::parse(&execution_policy) {
                Ok(parsed) => aot_config.execution_policy = parsed,
                Err(error) => report.push_warning(error),
            }
        }
        if let Some(value) = usize_field(section, "aot_residual_target_chunks") {
            aot_config.residual_target_chunks = Some(value);
        }
        if let Some(value) = usize_field(section, "aot_jacobian_target_chunks") {
            aot_config.jacobian_target_chunks = Some(value);
        }

        backend = backend.with_aot_config(aot_config);
    }

    backend
}

fn parse_postprocessing_config(
    postprocessing: Option<&HashMap<String, Option<Vec<Value>>>>,
    legacy_solver_section: Option<&HashMap<String, Option<Vec<Value>>>>,
    _report: &mut BvpGuiMigrationReport,
) -> BvpPostprocessingConfig {
    let mut config = BvpPostprocessingConfig::default();
    let section = postprocessing
        .filter(|section| !section.is_empty())
        .or(legacy_solver_section);

    if let Some(section) = section {
        if let Some(value) = bool_field(section, "plot") {
            config.plot = value;
        }
        if let Some(value) = bool_field(section, "gnuplot") {
            config.gnuplot = value;
        }
        if let Some(value) = bool_field(section, "save_to_csv") {
            config.save_to_csv = value;
        }
        if let Some(value) = string_field(section, "filename") {
            config.filename = value;
        }
        if let Some(value) = bool_field(section, "save") {
            config.save = value;
        }
        if let Some(value) = bool_field(section, "return_to_dimension") {
            config.return_to_dimension = value;
        }
        if let Some(value) = bool_field(section, "no_plots_in_terminal") {
            config.no_plots_in_terminal = value;
        }
        if let Some(value) = bool_field(section, "gui_plot") {
            config.gui_plot = value;
        }
    }
    config
}

fn parse_advanced_solver_config(
    rel_tolerance: Option<&HashMap<String, Option<Vec<Value>>>>,
    strategy_params: Option<&HashMap<String, Option<Vec<Value>>>>,
    report: &mut BvpGuiMigrationReport,
) -> BvpAdvancedSolverConfig {
    let mut config = BvpAdvancedSolverConfig::default();

    if let Some(section) = rel_tolerance {
        if let Some(value) = float_field(section, "C") {
            if value.is_finite() {
                config.rel_tolerance.c = value;
            } else {
                report.push_warning(
                    "Non-finite rel_tolerance.C encountered. Falling back to default.",
                );
            }
        }
        if let Some(value) = float_field(section, "J") {
            if value.is_finite() {
                config.rel_tolerance.j = value;
            } else {
                report.push_warning(
                    "Non-finite rel_tolerance.J encountered. Falling back to default.",
                );
            }
        }
        if let Some(value) = float_field(section, "Teta") {
            if value.is_finite() {
                config.rel_tolerance.teta = value;
            } else {
                report.push_warning(
                    "Non-finite rel_tolerance.Teta encountered. Falling back to default.",
                );
            }
        }
        if let Some(value) = float_field(section, "q") {
            if value.is_finite() {
                config.rel_tolerance.q = value;
            } else {
                report.push_warning(
                    "Non-finite rel_tolerance.q encountered. Falling back to default.",
                );
            }
        }
    }

    if let Some(section) = strategy_params {
        if let Some(value) = optional_usize_field(section, "max_jac") {
            config.strategy_params.max_jac = Some(value);
        }
        if let Some(value) = optional_usize_field(section, "max_damp_iter") {
            config.strategy_params.max_damp_iter = Some(value);
        }
        if let Some(value) = optional_float_field(section, "damp_factor") {
            if value.is_finite() {
                config.strategy_params.damp_factor = Some(value);
            } else {
                report.push_warning(
                    "Non-finite strategy_params.damp_factor encountered. Falling back to default.",
                );
            }
        }
    }

    config
}

fn parse_initial_guess_config(
    initial_guess: Option<&HashMap<String, Option<Vec<Value>>>>,
    report: &mut BvpGuiMigrationReport,
) -> BvpInitialGuessConfig {
    let mut config = BvpInitialGuessConfig::default();
    let Some(section) = initial_guess else {
        return config;
    };

    if let Some(value) = float_field(section, "universal") {
        if value.is_finite() {
            config.mode = BvpInitialGuessMode::Universal;
            config.universal = value;
            return config;
        }
        report.push_warning(
            "Non-finite initial_guess.universal encountered. Falling back to default.",
        );
    }

    let has_explicit_fields = section.contains_key("C")
        || section.contains_key("J")
        || section.contains_key("Teta")
        || section.contains_key("q");
    if has_explicit_fields {
        config.mode = BvpInitialGuessMode::PerVariable;
        if let Some(value) = float_field(section, "C") {
            config.c = value;
        } else {
            report.push_warning("Missing initial_guess.C. Falling back to default.");
        }
        if let Some(value) = float_field(section, "J") {
            config.j = value;
        } else {
            report.push_warning("Missing initial_guess.J. Falling back to default.");
        }
        if let Some(value) = float_field(section, "Teta") {
            config.teta = value;
        } else {
            report.push_warning("Missing initial_guess.Teta. Falling back to default.");
        }
        if let Some(value) = float_field(section, "q") {
            config.q = value;
        } else {
            report.push_warning("Missing initial_guess.q. Falling back to default.");
        }
    }

    config
}

fn parse_adaptive_grid_config(
    document: &DocumentMap,
    adaptive_strategy: Option<&HashMap<String, Option<Vec<Value>>>>,
    grid_refinement: Option<&HashMap<String, Option<Vec<Value>>>>,
    report: &mut BvpGuiMigrationReport,
) -> BvpAdaptiveGridConfig {
    let mut config = BvpAdaptiveGridConfig::default();
    let legacy_marker = document
        .get("strategy_params")
        .and_then(|section| section.get("adaptive"))
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(|value| match value {
            Value::Optional(None) => Some(false),
            Value::Optional(Some(inner)) => inner.as_boolean(),
            Value::String(value) if value.trim().eq_ignore_ascii_case("none") => Some(false),
            Value::String(value) if value.trim().eq_ignore_ascii_case("true") => Some(true),
            Value::String(value) if value.trim().eq_ignore_ascii_case("false") => Some(false),
            Value::Boolean(value) => Some(*value),
            _ => None,
        });

    let has_native_payload = adaptive_strategy.is_some_and(|section| !section.is_empty())
        || grid_refinement.is_some_and(|section| !section.is_empty());
    config.enabled = legacy_marker.unwrap_or(has_native_payload);

    if !config.enabled {
        return config;
    }

    if let Some(section) = adaptive_strategy {
        if let Some(value) = usize_field(section, "version") {
            if value != 1 {
                report.push_warning(format!(
                    "Adaptive strategy version `{value}` is not supported. Falling back to version 1."
                ));
            }
        }
        if let Some(value) = usize_field(section, "max_refinements") {
            config.max_refinements = value;
        }
    }

    if let Some(section) = grid_refinement {
        if let Some((method, params)) = parse_grid_refinement_section(section, report) {
            config.grid_refinement = BvpGridRefinementConfig { method, params };
        }
    }

    config
}

fn parse_grid_refinement_section(
    section: &HashMap<String, Option<Vec<Value>>>,
    report: &mut BvpGuiMigrationReport,
) -> Option<(BvpGridRefinementMethod, Vec<f64>)> {
    for method in [
        BvpGridRefinementMethod::Pearson,
        BvpGridRefinementMethod::GrcarSmooke,
        BvpGridRefinementMethod::Twopnt,
        BvpGridRefinementMethod::Easy,
        BvpGridRefinementMethod::DoublePoints,
    ] {
        if let Some((_raw_key, slot)) = section.iter().find(|(key, _)| {
            BvpGridRefinementMethod::parse(key).is_some_and(|parsed| parsed == method)
        }) {
            let params = slot
                .as_ref()
                .and_then(|values| {
                    if let Some(vector) = values.first().and_then(|value| value.as_vector()) {
                        Some(vector.clone())
                    } else {
                        let mut params = Vec::with_capacity(values.len());
                        for value in values {
                            params.push(value.as_float()?);
                        }
                        Some(params)
                    }
                })
                .unwrap_or_else(|| method.default_params());
            return Some((method, params));
        }
    }

    report.push_warning(
        "Unknown grid refinement method found. Falling back to `pearson`.".to_string(),
    );
    None
}

fn string_field(section: &HashMap<String, Option<Vec<Value>>>, field_name: &str) -> Option<String> {
    let value = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())?;

    match value {
        Value::String(value) => Some(value.clone()),
        Value::Optional(Some(inner)) => inner.as_string().cloned(),
        Value::Optional(None) => None,
        _ => None,
    }
}

fn optional_string_field(
    section: &HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) -> Option<String> {
    let value = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())?;

    match value {
        Value::Optional(None) => None,
        Value::Optional(Some(inner)) => inner.as_string().cloned(),
        Value::String(value) if value.eq_ignore_ascii_case("none") => None,
        Value::String(value) => Some(value.clone()),
        _ => None,
    }
}

fn bool_field(section: &HashMap<String, Option<Vec<Value>>>, field_name: &str) -> Option<bool> {
    let value = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())?;

    match value {
        Value::Boolean(value) => Some(*value),
        Value::String(value) if value.eq_ignore_ascii_case("true") => Some(true),
        Value::String(value) if value.eq_ignore_ascii_case("false") => Some(false),
        _ => None,
    }
}

fn optional_usize_field(
    section: &HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) -> Option<usize> {
    let value = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())?;

    match value {
        Value::Optional(None) => None,
        Value::Optional(Some(inner)) => try_value_as_usize(inner).ok(),
        Value::Integer(_) | Value::Usize(_) => try_value_as_usize(value).ok(),
        _ => None,
    }
}

fn optional_float_field(
    section: &HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) -> Option<f64> {
    let value = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())?;

    match value {
        Value::Optional(None) => None,
        Value::Optional(Some(inner)) => inner.as_float(),
        _ => value.as_float(),
    }
}

fn usize_field(section: &HashMap<String, Option<Vec<Value>>>, field_name: &str) -> Option<usize> {
    section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(|value| try_value_as_usize(value).ok())
}

fn float_field(section: &HashMap<String, Option<Vec<Value>>>, field_name: &str) -> Option<f64> {
    section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(Value::as_float)
}

fn set_string_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: &str,
) {
    section.insert(
        field_name.to_string(),
        Some(vec![Value::String(value.to_string())]),
    );
}

fn set_optional_string_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: Option<&str>,
) {
    let value = value
        .map(|value| Value::Optional(Some(Box::new(Value::String(value.to_string())))))
        .unwrap_or(Value::Optional(None));
    section.insert(field_name.to_string(), Some(vec![value]));
}

fn set_bool_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: bool,
) {
    section.insert(field_name.to_string(), Some(vec![Value::Boolean(value)]));
}

fn set_optional_usize_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: Option<usize>,
) {
    match value {
        Some(value) => {
            section.insert(
                field_name.to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::Usize(value))))]),
            );
        }
        None => {
            section.remove(field_name);
        }
    }
}

fn set_optional_float_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: Option<f64>,
) {
    match value {
        Some(value) => {
            section.insert(
                field_name.to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::Float(value))))]),
            );
        }
        None => {
            section.remove(field_name);
        }
    }
}

fn set_float_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: f64,
) {
    section.insert(field_name.to_string(), Some(vec![Value::Float(value)]));
}

fn set_usize_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    value: usize,
) {
    section.insert(field_name.to_string(), Some(vec![Value::Usize(value)]));
}

fn ensure_section_mut<'a>(
    document: &'a mut DocumentMap,
    section_name: &str,
) -> &'a mut HashMap<String, Option<Vec<Value>>> {
    document
        .entry(section_name.to_string())
        .or_insert_with(HashMap::new)
}

fn normalize_key(value: &str) -> String {
    value.trim().to_ascii_lowercase().replace('-', "_")
}

fn generated_backend_name(config: &ReactorBvpSolverConfig) -> String {
    let matrix_prefix = match config.matrix_backend {
        ReactorBvpMatrixBackend::Banded => "banded",
        ReactorBvpMatrixBackend::Sparse => "sparse",
    };

    match config.execution_backend {
        ReactorBvpExecutionBackend::Lambdify => format!("{}_lambdify", matrix_prefix),
        ReactorBvpExecutionBackend::Aot(aot) => match aot.compiler {
            ReactorBvpAotCompiler::CTcc => format!("{}_aot_tcc", matrix_prefix),
            ReactorBvpAotCompiler::CGcc => format!("{}_aot_gcc", matrix_prefix),
            ReactorBvpAotCompiler::Zig => format!("{}_aot_zig", matrix_prefix),
        },
    }
}

const BVP_AOT_ONLY_FIELDS: &[&str] = &[
    "aot_codegen_backend",
    "aot_c_compiler",
    "aot_build_policy",
    "aot_build_profile",
    "aot_compile_preset",
    "aot_execution_policy",
    "aot_residual_target_chunks",
    "aot_jacobian_target_chunks",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum GuiExecutionBackend {
    Lambdify,
    Aot,
}

impl GuiExecutionBackend {
    pub(crate) fn from_generated_backend(raw: &str) -> Option<Self> {
        let normalized = raw.trim().to_ascii_lowercase();
        if normalized.contains("lambdify") {
            Some(Self::Lambdify)
        } else if generated_backend_requests_aot(&normalized) {
            Some(Self::Aot)
        } else {
            None
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum GuiAotCodegenBackend {
    C,
    Zig,
    Rust,
}

impl GuiAotCodegenBackend {
    pub(crate) fn parse(raw: &str) -> Option<Self> {
        match normalize_key(raw).as_str() {
            "c" => Some(Self::C),
            "zig" => Some(Self::Zig),
            "rust" | "rs" => Some(Self::Rust),
            _ => None,
        }
    }

    pub(crate) fn infer_from_preset(preset: &str) -> Self {
        let normalized = preset.trim().to_ascii_lowercase();
        if normalized.contains("zig") {
            Self::Zig
        } else if normalized.contains("rust") {
            Self::Rust
        } else {
            Self::C
        }
    }

    pub(crate) fn as_document_value(self) -> &'static str {
        match self {
            Self::C => "C",
            Self::Zig => "Zig",
            Self::Rust => "Rust",
        }
    }
}

fn ensure_field_slot<'a>(
    section: &'a mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    default_value: Value,
) -> &'a mut Option<Vec<Value>> {
    let entry = section
        .entry(field_name.to_string())
        .or_insert_with(|| Some(vec![default_value.clone()]));

    if entry.is_none() {
        *entry = Some(vec![default_value.clone()]);
    }

    if entry.as_ref().is_some_and(|values| values.is_empty()) {
        *entry = Some(vec![default_value]);
    }

    entry
}

pub(crate) fn generated_backend_requests_aot(generated_backend: &str) -> bool {
    let normalized = generated_backend.trim().to_ascii_lowercase();
    normalized.contains("_aot")
        || normalized.starts_with("aot_")
        || normalized.ends_with("_build_if_missing")
        || normalized.ends_with("_require_prebuilt")
        || normalized.ends_with("_rebuild_always")
        || normalized.ends_with("_repeated")
        || normalized.contains("_atomview_gcc")
        || normalized.contains("_atomview_tcc")
        || normalized.contains("_atomview_zig")
}

fn compiler_inferred_from_preset(preset: &str) -> &'static str {
    let normalized = preset.trim().to_ascii_lowercase();
    if normalized.ends_with("_gcc") {
        "gcc"
    } else if normalized.ends_with("_zig") {
        "zig"
    } else {
        "tcc"
    }
}

/// Normalizes the visible solver-backend controls into a single canonical raw section.
///
/// The GUI and the reactor-task parser share this helper so the renderer does not
/// own a second copy of the backend contract. When a document only carries the
/// compatibility `backend_policy`, the helper still reconstructs the visible
/// backend mode instead of forcing callers to duplicate that fallback.
pub(crate) fn normalize_bvp_solver_backend_section(
    section: &mut HashMap<String, Option<Vec<Value>>>,
) {
    let generated_backend = string_field(section, "generated_backend");
    let backend_mode = if let Some(raw_backend) = generated_backend.as_deref() {
        match GuiExecutionBackend::from_generated_backend(raw_backend) {
            Some(mode) => mode,
            None => return,
        }
    } else {
        match string_field(section, "backend_policy")
            .map(|value| normalize_key(&value))
            .as_deref()
        {
            Some("aot_only") => GuiExecutionBackend::Aot,
            Some("lambdify_only") => GuiExecutionBackend::Lambdify,
            _ => return,
        }
    };

    let matrix_backend = string_field(section, "matrix_backend")
        .or_else(|| string_field(section, "method"))
        .unwrap_or_else(|| "Banded".to_string());
    let matrix_prefix = if matrix_backend.trim().eq_ignore_ascii_case("Sparse") {
        "sparse"
    } else {
        "banded"
    };

    match backend_mode {
        GuiExecutionBackend::Lambdify => {
            set_string_field(section, "backend_policy", "lambdify_only");
            for field_name in BVP_AOT_ONLY_FIELDS {
                section.remove(*field_name);
            }

            let next_backend = format!("{}_lambdify", matrix_prefix);
            if generated_backend.as_deref() != Some(next_backend.as_str()) {
                set_string_field(section, "generated_backend", &next_backend);
                info!("Field 'generated_backend' changed to: {}", next_backend);
            }
        }
        GuiExecutionBackend::Aot => {
            set_string_field(section, "backend_policy", "aot_only");

            let codegen = string_field(section, "aot_codegen_backend")
                .as_deref()
                .and_then(GuiAotCodegenBackend::parse)
                .unwrap_or_else(|| {
                    GuiAotCodegenBackend::infer_from_preset(
                        generated_backend.as_deref().unwrap_or_default(),
                    )
                });
            set_string_field(section, "aot_codegen_backend", codegen.as_document_value());

            match codegen {
                GuiAotCodegenBackend::C => {
                    if !section.contains_key("aot_c_compiler") {
                        set_string_field(
                            section,
                            "aot_c_compiler",
                            compiler_inferred_from_preset(
                                generated_backend.as_deref().unwrap_or_default(),
                            ),
                        );
                    }
                }
                _ => {
                    section.remove("aot_c_compiler");
                }
            }

            for (field_name, default_value) in [
                (
                    "aot_build_policy",
                    Value::String("build_if_missing".to_string()),
                ),
                ("aot_build_profile", Value::String("release".to_string())),
                (
                    "aot_compile_preset",
                    Value::String("dev_fastest".to_string()),
                ),
                ("aot_execution_policy", Value::String("auto".to_string())),
            ] {
                let _ = ensure_field_slot(section, field_name, default_value);
            }

            let next_backend = match codegen {
                GuiAotCodegenBackend::C => {
                    let compiler = string_field(section, "aot_c_compiler").unwrap_or_else(|| {
                        compiler_inferred_from_preset(
                            generated_backend.as_deref().unwrap_or_default(),
                        )
                        .to_string()
                    });
                    format!("{}_aot_{}", matrix_prefix, compiler.to_ascii_lowercase())
                }
                GuiAotCodegenBackend::Zig => format!("{}_aot_zig", matrix_prefix),
                GuiAotCodegenBackend::Rust => format!("{}_aot", matrix_prefix),
            };
            if generated_backend.as_deref() != Some(next_backend.as_str()) {
                set_string_field(section, "generated_backend", &next_backend);
                info!("Field 'generated_backend' changed to: {}", next_backend);
            }
        }
    }
}

/// Executables that the GUI permits for AOT compilation.
///
/// RustedSciThe accepts explicit compiler names, but the GUI should keep a
/// conservative allow-list so shared documents cannot smuggle an arbitrary
/// executable into an AOT confirmation flow.
const BVP_APPROVED_AOT_COMPILERS: &[&str] = &["tcc", "gcc"];

/// Stable description of an AOT toolchain request extracted from solver settings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AotToolchainRequest {
    pub generated_backend: String,
    pub codegen_backend: String,
    pub compiler: Option<String>,
    pub build_policy: String,
}

impl AotToolchainRequest {
    pub(crate) fn summary(&self) -> String {
        format!(
            "backend={}, codegen={}, compiler={}, build_policy={}",
            self.generated_backend,
            self.codegen_backend,
            self.compiler.as_deref().unwrap_or("backend default"),
            self.build_policy
        )
    }
}

/// Validate the GUI-visible AOT toolchain fields before solver handoff.
pub(crate) fn validate_aot_toolchain_fields(document: &DocumentMap) -> Result<(), String> {
    let Some(section) = document.get("solver_settings") else {
        return Ok(());
    };

    let codegen = string_field(section, "aot_codegen_backend")
        .map(|raw| {
            GuiAotCodegenBackend::parse(&raw).ok_or_else(|| {
                format!(
                    "AOT codegen backend `{}` is invalid. Choose C, Zig, or Rust",
                    raw
                )
            })
        })
        .transpose()?;
    if let Some(compiler) = string_field(section, "aot_c_compiler") {
        let normalized = compiler.trim().to_ascii_lowercase();
        if !BVP_APPROVED_AOT_COMPILERS.contains(&normalized.as_str()) {
            return Err(format!(
                "AOT compiler `{}` is not allowed by the GUI. Choose one of: {}",
                compiler,
                BVP_APPROVED_AOT_COMPILERS.join(", ")
            ));
        }
        if let Some(non_c_codegen) = codegen.filter(|codegen| *codegen != GuiAotCodegenBackend::C) {
            return Err(format!(
                "AOT C compiler `{}` cannot be used with the {} codegen backend",
                compiler,
                non_c_codegen.as_document_value()
            ));
        }
    }

    Ok(())
}

/// Build the AOT confirmation payload only when the effective document requests it.
pub(crate) fn aot_toolchain_request(document: &DocumentMap) -> Option<AotToolchainRequest> {
    let section = document.get("solver_settings")?;
    let generated_backend = string_field(section, "generated_backend")?;
    let policy_requests_aot = string_field(section, "backend_policy")
        .is_some_and(|policy| policy.trim().eq_ignore_ascii_case("aot_only"));
    if !generated_backend_requests_aot(&generated_backend) && !policy_requests_aot {
        return None;
    }

    Some(AotToolchainRequest {
        generated_backend,
        codegen_backend: string_field(section, "aot_codegen_backend")
            .unwrap_or_else(|| "C".to_string()),
        compiler: string_field(section, "aot_c_compiler"),
        build_policy: string_field(section, "aot_build_policy")
            .unwrap_or_else(|| "build_if_missing".to_string()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ReactorsBVP::task_parser_reactor_BVP::SIMPLE_BVP_TEMPLATE;
    use RustedSciThe::command_interpreter::task_parser::DocumentParser;

    fn parse_template() -> DocumentMap {
        let mut parser = DocumentParser::new(SIMPLE_BVP_TEMPLATE.to_string());
        parser.parse_document().expect("template must parse");
        parser
            .get_result()
            .cloned()
            .expect("template document must exist")
    }

    fn section_value<'a>(document: &'a DocumentMap, section: &str, field: &str) -> &'a Value {
        document
            .get(section)
            .and_then(|section| section.get(field))
            .and_then(|slot| slot.as_ref())
            .and_then(|values| values.first())
            .expect("expected field to exist")
    }

    #[test]
    fn default_gui_config_roundtrips_canonical_solver_settings() {
        let document = parse_template();
        let (config, report) = BvpGuiConfig::from_document(&document);

        assert!(
            report.warnings.is_empty(),
            "warnings: {:?}",
            report.warnings
        );
        assert_eq!(config.core.scheme, BvpScheme::Forward);
        assert_eq!(
            config.backend.matrix_backend,
            ReactorBvpMatrixBackend::Banded
        );
        assert_eq!(
            config.backend.symbolic_backend,
            ReactorBvpSymbolicBackend::AtomView
        );
        assert!(matches!(
            config.backend.execution_backend,
            ReactorBvpExecutionBackend::Lambdify
        ));
        assert!(config.adaptive_grid.enabled);
        assert_eq!(config.adaptive_grid.version, 1);
        assert_eq!(
            config.adaptive_grid.grid_refinement.method,
            BvpGridRefinementMethod::GrcarSmooke
        );
        assert_eq!(config.initial_guess.mode, BvpInitialGuessMode::Universal);
        assert_eq!(config.initial_guess.universal, 1e-2);
        assert_eq!(config.advanced_solver.rel_tolerance.c, 1e-7);
        assert_eq!(config.advanced_solver.strategy_params.max_jac, Some(3));

        let mut rebuilt = DocumentMap::new();
        config.apply_to_document(&mut rebuilt);

        assert_eq!(
            section_value(&rebuilt, "solver_settings", "scheme"),
            &Value::String("forward".to_string())
        );
        assert_eq!(
            section_value(&rebuilt, "solver_settings", "matrix_backend"),
            &Value::String("Banded".to_string())
        );
        assert_eq!(
            section_value(&rebuilt, "solver_settings", "generated_backend"),
            &Value::String("banded_lambdify".to_string())
        );
        assert_eq!(
            section_value(&rebuilt, "initial_guess", "universal"),
            &Value::Float(1e-2)
        );
        assert_eq!(
            section_value(&rebuilt, "rel_tolerance", "C"),
            &Value::Float(1e-7)
        );
    }

    #[test]
    fn aot_backend_roundtrips_through_typed_config() {
        let mut document = parse_template();
        let solver = document
            .get_mut("solver_settings")
            .expect("solver_settings");
        set_string_field(solver, "generated_backend", "sparse_aot_gcc");
        set_string_field(solver, "matrix_backend", "Sparse");
        set_string_field(solver, "symbolic_backend", "ExprLegacy");
        set_string_field(solver, "aot_codegen_backend", "C");
        set_string_field(solver, "aot_c_compiler", "gcc");
        set_string_field(solver, "aot_build_policy", "rebuild_always");
        set_string_field(solver, "aot_build_profile", "debug");
        set_string_field(solver, "aot_compile_preset", "fast_build");
        set_string_field(solver, "aot_execution_policy", "parallel");
        set_usize_field(solver, "aot_residual_target_chunks", 12);
        set_usize_field(solver, "aot_jacobian_target_chunks", 9);

        let (config, report) = BvpGuiConfig::from_document(&document);
        assert!(
            report.warnings.is_empty(),
            "warnings: {:?}",
            report.warnings
        );
        assert_eq!(
            config.backend.matrix_backend,
            ReactorBvpMatrixBackend::Sparse
        );
        assert_eq!(
            config.backend.symbolic_backend,
            ReactorBvpSymbolicBackend::ExprLegacy
        );
        match config.backend.execution_backend {
            ReactorBvpExecutionBackend::Aot(aot) => {
                assert_eq!(aot.compiler, ReactorBvpAotCompiler::CGcc);
                assert_eq!(aot.build_policy, ReactorBvpAotBuildPolicy::RebuildAlways);
                assert_eq!(aot.build_profile, ReactorBvpAotBuildProfile::Debug);
                assert_eq!(
                    aot.compile_preset,
                    crate::ReactorsBVP::solver_backend::ReactorBvpAotCompilePreset::FastBuild
                );
                assert_eq!(aot.execution_policy, ReactorBvpAotExecutionPolicy::Parallel);
                assert_eq!(aot.residual_target_chunks, Some(12));
                assert_eq!(aot.jacobian_target_chunks, Some(9));
            }
            other => panic!("expected AOT backend, got {other:?}"),
        }

        let mut rebuilt = DocumentMap::new();
        config.apply_to_document(&mut rebuilt);
        assert_eq!(
            section_value(&rebuilt, "solver_settings", "generated_backend"),
            &Value::String("sparse_aot_gcc".to_string())
        );
        assert_eq!(
            section_value(&rebuilt, "solver_settings", "backend_policy"),
            &Value::String("aot_only".to_string())
        );
    }

    #[test]
    fn aot_toolchain_request_is_extracted_from_shared_backend_state() {
        let mut document = parse_template();
        let solver = document
            .get_mut("solver_settings")
            .expect("solver_settings");
        set_string_field(solver, "generated_backend", "banded_aot_gcc");
        set_string_field(solver, "aot_codegen_backend", "C");
        set_string_field(solver, "aot_c_compiler", "gcc");
        set_string_field(solver, "aot_build_policy", "build_if_missing");

        let request = aot_toolchain_request(&document).expect("expected AOT request");
        assert_eq!(request.generated_backend, "banded_aot_gcc");
        assert_eq!(request.codegen_backend, "C");
        assert_eq!(request.compiler.as_deref(), Some("gcc"));
        assert_eq!(request.build_policy, "build_if_missing");
        assert!(request.summary().contains("backend=banded_aot_gcc"));
    }

    #[test]
    fn initial_guess_roundtrips_through_typed_config() {
        let mut document = parse_template();
        document.insert(
            "initial_guess".to_string(),
            HashMap::from([
                ("C".to_string(), Some(vec![Value::Float(0.2)])),
                ("J".to_string(), Some(vec![Value::Float(0.3)])),
                ("Teta".to_string(), Some(vec![Value::Float(0.4)])),
                ("q".to_string(), Some(vec![Value::Float(0.5)])),
            ]),
        );

        let (config, report) = BvpGuiConfig::from_document(&document);
        assert!(
            report.warnings.is_empty(),
            "warnings: {:?}",
            report.warnings
        );
        assert!(matches!(
            config.initial_guess.mode,
            BvpInitialGuessMode::PerVariable
        ));
        assert_eq!(config.initial_guess.c, 0.2);
        assert_eq!(config.initial_guess.j, 0.3);
        assert_eq!(config.initial_guess.teta, 0.4);
        assert_eq!(config.initial_guess.q, 0.5);

        let mut rebuilt = DocumentMap::new();
        config.apply_to_document(&mut rebuilt);
        assert_eq!(
            section_value(&rebuilt, "initial_guess", "C"),
            &Value::Float(0.2)
        );
        assert_eq!(
            section_value(&rebuilt, "initial_guess", "q"),
            &Value::Float(0.5)
        );
        assert!(
            !rebuilt
                .get("initial_guess")
                .expect("initial_guess section")
                .contains_key("universal")
        );
    }

    #[test]
    fn advanced_solver_roundtrips_through_typed_config() {
        let mut document = parse_template();
        document.insert(
            "rel_tolerance".to_string(),
            HashMap::from([
                ("C".to_string(), Some(vec![Value::Float(1e-5)])),
                ("J".to_string(), Some(vec![Value::Float(2e-5)])),
                ("Teta".to_string(), Some(vec![Value::Float(3e-5)])),
                ("q".to_string(), Some(vec![Value::Float(4e-5)])),
            ]),
        );
        document.insert(
            "strategy_params".to_string(),
            HashMap::from([
                (
                    "max_jac".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Usize(7))))]),
                ),
                (
                    "max_damp_iter".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Usize(11))))]),
                ),
                (
                    "damp_factor".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Float(0.25))))]),
                ),
            ]),
        );

        let (config, report) = BvpGuiConfig::from_document(&document);
        assert!(
            report.warnings.is_empty(),
            "warnings: {:?}",
            report.warnings
        );
        assert_eq!(config.advanced_solver.rel_tolerance.c, 1e-5);
        assert_eq!(config.advanced_solver.rel_tolerance.q, 4e-5);
        assert_eq!(config.advanced_solver.strategy_params.max_jac, Some(7));
        assert_eq!(
            config.advanced_solver.strategy_params.max_damp_iter,
            Some(11)
        );
        assert_eq!(
            config.advanced_solver.strategy_params.damp_factor,
            Some(0.25)
        );

        let mut rebuilt = DocumentMap::new();
        config.apply_to_document(&mut rebuilt);
        assert_eq!(
            section_value(&rebuilt, "rel_tolerance", "J"),
            &Value::Float(2e-5)
        );
        assert_eq!(
            section_value(&rebuilt, "strategy_params", "max_damp_iter"),
            &Value::Optional(Some(Box::new(Value::Usize(11))))
        );
    }

    #[test]
    fn validate_aot_toolchain_fields_rejects_unapproved_compiler() {
        let mut document = parse_template();
        let solver = document
            .get_mut("solver_settings")
            .expect("solver_settings");
        set_string_field(solver, "generated_backend", "banded_aot_tcc");
        set_string_field(solver, "aot_codegen_backend", "C");
        set_string_field(solver, "aot_c_compiler", "clang");

        let error = validate_aot_toolchain_fields(&document)
            .expect_err("unapproved compiler must be rejected");
        assert!(error.contains("not allowed"));
    }

    #[test]
    fn legacy_adaptive_marker_becomes_typed_adaptive_state() {
        let mut document = parse_template();
        let strategy_params = document
            .get_mut("strategy_params")
            .expect("strategy_params section");
        strategy_params.insert("adaptive".to_string(), Some(vec![Value::Optional(None)]));

        let (config, report) = BvpGuiConfig::from_document(&document);
        assert!(
            report.warnings.is_empty(),
            "warnings: {:?}",
            report.warnings
        );
        assert!(!config.adaptive_grid.enabled);

        let mut rebuilt = document.clone();
        config.apply_to_document(&mut rebuilt);
        assert!(!rebuilt.contains_key("adaptive_strategy"));
        assert!(!rebuilt.contains_key("grid_refinement"));
    }

    #[test]
    fn typed_snapshot_apply_is_idempotent_on_canonical_documents() {
        let document = parse_template();
        let (config, report) = BvpGuiConfig::from_document(&document);
        assert!(report.warnings.is_empty());

        let mut once = document.clone();
        config.apply_to_document(&mut once);

        let mut twice = once.clone();
        config.apply_to_document(&mut twice);

        assert_eq!(once, twice);
    }

    #[test]
    fn normalize_backend_section_keeps_aot_codegen_choice_and_updates_backend_name() {
        let mut document = parse_template();
        {
            let solver = document
                .get_mut("solver_settings")
                .expect("solver_settings");
            set_string_field(solver, "generated_backend", "banded_aot");
            set_string_field(solver, "matrix_backend", "Banded");
            set_string_field(solver, "aot_codegen_backend", "Rust");
            set_string_field(solver, "backend_policy", "lambdify_only");

            normalize_bvp_solver_backend_section(solver);

            assert!(!solver.contains_key("aot_c_compiler"));

            set_string_field(solver, "matrix_backend", "Sparse");
            normalize_bvp_solver_backend_section(solver);
        }

        assert_eq!(
            section_value(&document, "solver_settings", "generated_backend"),
            &Value::String("sparse_aot".to_string())
        );
        assert_eq!(
            section_value(&document, "solver_settings", "backend_policy"),
            &Value::String("aot_only".to_string())
        );
        assert_eq!(
            section_value(&document, "solver_settings", "aot_codegen_backend"),
            &Value::String("Rust".to_string())
        );
    }

    #[test]
    fn normalize_backend_section_removes_aot_fields_for_lambdify() {
        let mut document = parse_template();
        {
            let solver = document
                .get_mut("solver_settings")
                .expect("solver_settings");
            set_string_field(solver, "generated_backend", "banded_lambdify");
            set_string_field(solver, "matrix_backend", "Banded");
            set_string_field(solver, "aot_codegen_backend", "C");
            set_string_field(solver, "aot_c_compiler", "gcc");
            set_string_field(solver, "aot_build_policy", "rebuild_always");
            set_string_field(solver, "aot_build_profile", "debug");

            normalize_bvp_solver_backend_section(solver);
        }

        assert_eq!(
            section_value(&document, "solver_settings", "generated_backend"),
            &Value::String("banded_lambdify".to_string())
        );
        assert_eq!(
            section_value(&document, "solver_settings", "backend_policy"),
            &Value::String("lambdify_only".to_string())
        );
        let solver = document.get("solver_settings").expect("solver_settings");
        assert!(!solver.contains_key("aot_codegen_backend"));
        assert!(!solver.contains_key("aot_c_compiler"));
        assert!(!solver.contains_key("aot_build_policy"));
        assert!(!solver.contains_key("aot_build_profile"));
    }

    #[test]
    fn normalize_backend_section_recovers_backend_policy_only_aot_docs() {
        let mut document = parse_template();
        {
            let solver = document
                .get_mut("solver_settings")
                .expect("solver_settings");
            solver.remove("generated_backend");
            set_string_field(solver, "method", "Sparse");
            set_string_field(solver, "matrix_backend", "Sparse");
            set_string_field(solver, "backend_policy", "aot_only");
            set_string_field(solver, "aot_c_compiler", "gcc");

            normalize_bvp_solver_backend_section(solver);
        }

        assert_eq!(
            section_value(&document, "solver_settings", "generated_backend"),
            &Value::String("sparse_aot_gcc".to_string())
        );
        assert_eq!(
            section_value(&document, "solver_settings", "backend_policy"),
            &Value::String("aot_only".to_string())
        );
        assert_eq!(
            section_value(&document, "solver_settings", "aot_c_compiler"),
            &Value::String("gcc".to_string())
        );
    }
}
