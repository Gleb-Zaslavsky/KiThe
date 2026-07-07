//! Solver-backend facade for reactor BVP tasks.
//!
//! KiThe reactor problems always start from symbolic equations, so this facade
//! intentionally exposes only the meaningful matrix/execution routes for that
//! workflow: sparse or banded generated callbacks, with lambdify as the
//! toolchain-free default and AOT available through explicit user settings.

use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::DampedSolverOptions;
use RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::{
    AotBuildPolicy as RustedAotBuildPolicy, AotBuildProfile as RustedAotBuildProfile,
    AotChunkingPolicy, AotExecutionPolicy as RustedAotExecutionPolicy, GeneratedBackendConfig,
};
use RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy;
use RustedSciThe::symbolic::codegen::codegen_orchestrator::ParallelExecutorConfig;
use RustedSciThe::symbolic::codegen::codegen_runtime_api::ResidualChunkingStrategy;
use RustedSciThe::symbolic::codegen::codegen_tasks::SparseChunkingStrategy;
use RustedSciThe::symbolic::symbolic_functions_BVP::BvpSymbolicAssemblyBackend;

/// Matrix representation used by generated reactor BVP callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpMatrixBackend {
    /// Native banded Jacobian representation. This is the KiThe default.
    Banded,
    /// Sparse Jacobian representation for compatibility and diagnostics.
    Sparse,
}

impl Default for ReactorBvpMatrixBackend {
    fn default() -> Self {
        Self::Banded
    }
}

impl ReactorBvpMatrixBackend {
    /// Parse user-facing matrix backend names used in reactor task documents.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "banded" => Ok(Self::Banded),
            "sparse" => Ok(Self::Sparse),
            other => Err(format!(
                "Unknown reactor BVP matrix backend `{other}`. Expected `banded` or `sparse`"
            )),
        }
    }
}

/// Build profile used by AOT artifact policies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpAotBuildProfile {
    Release,
    Debug,
}

impl Default for ReactorBvpAotBuildProfile {
    fn default() -> Self {
        Self::Release
    }
}

impl ReactorBvpAotBuildProfile {
    /// Parse user-facing AOT build profile names.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "release" => Ok(Self::Release),
            "debug" => Ok(Self::Debug),
            other => Err(format!(
                "Unknown reactor BVP AOT build profile `{other}`. Expected `release` or `debug`"
            )),
        }
    }

    fn to_rusted(self) -> RustedAotBuildProfile {
        match self {
            Self::Release => RustedAotBuildProfile::Release,
            Self::Debug => RustedAotBuildProfile::Debug,
        }
    }
}

/// Build lifecycle policy for generated AOT artifacts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpAotBuildPolicy {
    UseIfAvailable,
    BuildIfMissing,
    RequirePrebuilt,
    RebuildAlways,
}

impl Default for ReactorBvpAotBuildPolicy {
    fn default() -> Self {
        Self::BuildIfMissing
    }
}

impl ReactorBvpAotBuildPolicy {
    /// Parse user-facing AOT build policy names.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "use_if_available" | "useavailable" => Ok(Self::UseIfAvailable),
            "build_if_missing" | "buildifmissing" => Ok(Self::BuildIfMissing),
            "require_prebuilt" | "requireprebuilt" => Ok(Self::RequirePrebuilt),
            "rebuild" | "rebuild_always" | "rebuildalways" => Ok(Self::RebuildAlways),
            other => Err(format!(
                "Unknown reactor BVP AOT build policy `{other}`. Expected `use_if_available`, `build_if_missing`, `require_prebuilt`, or `rebuild_always`"
            )),
        }
    }

    fn to_rusted(self, profile: ReactorBvpAotBuildProfile) -> RustedAotBuildPolicy {
        match self {
            Self::UseIfAvailable => RustedAotBuildPolicy::UseIfAvailable,
            Self::BuildIfMissing => RustedAotBuildPolicy::BuildIfMissing {
                profile: profile.to_rusted(),
            },
            Self::RequirePrebuilt => RustedAotBuildPolicy::RequirePrebuilt,
            Self::RebuildAlways => RustedAotBuildPolicy::RebuildAlways {
                profile: profile.to_rusted(),
            },
        }
    }
}

/// Compile-time preset for generated AOT artifacts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpAotCompilePreset {
    Production,
    FastBuild,
    DevFastest,
}

impl Default for ReactorBvpAotCompilePreset {
    fn default() -> Self {
        Self::DevFastest
    }
}

impl ReactorBvpAotCompilePreset {
    /// Parse user-facing AOT compile preset names.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "production" | "prod" => Ok(Self::Production),
            "fast_build" | "fastbuild" => Ok(Self::FastBuild),
            "dev_fastest" | "devfastest" | "fastest" => Ok(Self::DevFastest),
            other => Err(format!(
                "Unknown reactor BVP AOT compile preset `{other}`. Expected `production`, `fast_build`, or `dev_fastest`"
            )),
        }
    }
}

/// Runtime execution policy for compiled AOT callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpAotExecutionPolicy {
    Auto,
    Sequential,
    Parallel,
}

impl Default for ReactorBvpAotExecutionPolicy {
    fn default() -> Self {
        Self::Auto
    }
}

impl ReactorBvpAotExecutionPolicy {
    /// Parse user-facing AOT execution policy names.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "auto" => Ok(Self::Auto),
            "sequential" | "sequential_only" => Ok(Self::Sequential),
            "parallel" => Ok(Self::Parallel),
            other => Err(format!(
                "Unknown reactor BVP AOT execution policy `{other}`. Expected `auto`, `sequential`, or `parallel`"
            )),
        }
    }

    fn to_rusted(self) -> RustedAotExecutionPolicy {
        match self {
            Self::Auto => RustedAotExecutionPolicy::Auto,
            Self::Sequential => RustedAotExecutionPolicy::SequentialOnly,
            Self::Parallel => RustedAotExecutionPolicy::Parallel(ParallelExecutorConfig::default()),
        }
    }
}

/// Symbolic assembly route used before lambdify/AOT lowering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpSymbolicBackend {
    /// Modern AtomView symbolic assembly path. This is the KiThe default.
    AtomView,
    /// Historical expression-tree assembly path for compatibility checks.
    ExprLegacy,
}

impl Default for ReactorBvpSymbolicBackend {
    fn default() -> Self {
        Self::AtomView
    }
}

impl ReactorBvpSymbolicBackend {
    /// Parse user-facing symbolic backend names used in reactor task documents.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "atomview" | "atom_view" => Ok(Self::AtomView),
            "exprlegacy" | "expr_legacy" | "legacy" => Ok(Self::ExprLegacy),
            other => Err(format!(
                "Unknown reactor BVP symbolic backend `{other}`. Expected `AtomView` or `ExprLegacy`"
            )),
        }
    }

    fn to_rusted(self) -> BvpSymbolicAssemblyBackend {
        match self {
            Self::AtomView => BvpSymbolicAssemblyBackend::AtomView,
            Self::ExprLegacy => BvpSymbolicAssemblyBackend::ExprLegacy,
        }
    }
}

/// AOT compiler choice for generated callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpAotCompiler {
    /// C backend compiled with TinyCC.
    CTcc,
    /// C backend compiled with GCC.
    CGcc,
    /// Zig backend.
    Zig,
}

impl Default for ReactorBvpAotCompiler {
    fn default() -> Self {
        Self::CTcc
    }
}

impl ReactorBvpAotCompiler {
    /// Parse user-facing AOT compiler names.
    pub fn parse(value: &str) -> Result<Self, String> {
        match normalize_key(value).as_str() {
            "tcc" | "ctcc" | "c_tcc" => Ok(Self::CTcc),
            "gcc" | "cgcc" | "c_gcc" => Ok(Self::CGcc),
            "zig" => Ok(Self::Zig),
            other => Err(format!(
                "Unknown reactor BVP AOT compiler `{other}`. Expected `tcc`, `gcc`, or `zig`"
            )),
        }
    }
}

/// AOT configuration for generated residual/Jacobian callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactorBvpAotConfig {
    /// Toolchain/backend used to produce compiled callbacks.
    pub compiler: ReactorBvpAotCompiler,
    /// Artifact lifecycle policy used before a solve tries to load callbacks.
    pub build_policy: ReactorBvpAotBuildPolicy,
    /// Release/debug profile used by build and rebuild policies.
    pub build_profile: ReactorBvpAotBuildProfile,
    /// Compilation preset forwarded to RustedSciThe codegen.
    pub compile_preset: ReactorBvpAotCompilePreset,
    /// Runtime scheduling policy for compiled callbacks.
    pub execution_policy: ReactorBvpAotExecutionPolicy,
    /// Optional residual chunk target for AOT lowering.
    pub residual_target_chunks: Option<usize>,
    /// Optional sparse-Jacobian chunk target for AOT lowering.
    pub jacobian_target_chunks: Option<usize>,
}

impl Default for ReactorBvpAotConfig {
    fn default() -> Self {
        Self {
            compiler: ReactorBvpAotCompiler::default(),
            build_policy: ReactorBvpAotBuildPolicy::default(),
            build_profile: ReactorBvpAotBuildProfile::default(),
            compile_preset: ReactorBvpAotCompilePreset::default(),
            execution_policy: ReactorBvpAotExecutionPolicy::default(),
            residual_target_chunks: None,
            jacobian_target_chunks: None,
        }
    }
}

impl ReactorBvpAotConfig {
    /// Build an AOT config with a selected compiler and default lifecycle policy.
    pub fn for_compiler(compiler: ReactorBvpAotCompiler) -> Self {
        Self {
            compiler,
            ..Self::default()
        }
    }

    /// Override the AOT compiler.
    pub fn with_compiler(mut self, compiler: ReactorBvpAotCompiler) -> Self {
        self.compiler = compiler;
        self
    }

    /// Override the AOT build policy.
    pub fn with_build_policy(mut self, build_policy: ReactorBvpAotBuildPolicy) -> Self {
        self.build_policy = build_policy;
        self
    }

    /// Override the AOT build profile used by build/rebuild policies.
    pub fn with_build_profile(mut self, build_profile: ReactorBvpAotBuildProfile) -> Self {
        self.build_profile = build_profile;
        self
    }

    /// Override the AOT compile preset.
    pub fn with_compile_preset(mut self, compile_preset: ReactorBvpAotCompilePreset) -> Self {
        self.compile_preset = compile_preset;
        self
    }

    /// Override runtime execution policy for compiled callbacks.
    pub fn with_execution_policy(mut self, execution_policy: ReactorBvpAotExecutionPolicy) -> Self {
        self.execution_policy = execution_policy;
        self
    }

    /// Override target chunk counts used during generated residual/Jacobian lowering.
    pub fn with_target_chunks(
        mut self,
        residual_target_chunks: Option<usize>,
        jacobian_target_chunks: Option<usize>,
    ) -> Self {
        self.residual_target_chunks = residual_target_chunks;
        self.jacobian_target_chunks = jacobian_target_chunks;
        self
    }
}

/// Execution backend for generated residual/Jacobian callbacks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorBvpExecutionBackend {
    /// In-process lambdified callbacks.
    Lambdify,
    /// Ahead-of-time compiled callbacks.
    Aot(ReactorBvpAotConfig),
}

impl Default for ReactorBvpExecutionBackend {
    fn default() -> Self {
        Self::Lambdify
    }
}

/// User-facing backend configuration for KiThe reactor BVP solves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReactorBvpSolverConfig {
    pub matrix_backend: ReactorBvpMatrixBackend,
    pub symbolic_backend: ReactorBvpSymbolicBackend,
    pub execution_backend: ReactorBvpExecutionBackend,
}

impl ReactorBvpSolverConfig {
    /// KiThe default: lambdify execution, AtomView symbolic assembly, banded matrix callbacks.
    pub fn default_lambdify() -> Self {
        Self::default()
    }

    /// Compatibility route for comparing against sparse generated callbacks.
    pub fn sparse_lambdify() -> Self {
        Self {
            matrix_backend: ReactorBvpMatrixBackend::Sparse,
            ..Self::default_lambdify()
        }
    }

    /// Banded AtomView AOT route with a selected compiler.
    pub fn banded_aot(compiler: ReactorBvpAotCompiler) -> Self {
        Self {
            matrix_backend: ReactorBvpMatrixBackend::Banded,
            execution_backend: ReactorBvpExecutionBackend::Aot(ReactorBvpAotConfig::for_compiler(
                compiler,
            )),
            ..Self::default_lambdify()
        }
    }

    /// Sparse AtomView AOT route with a selected compiler.
    pub fn sparse_aot(compiler: ReactorBvpAotCompiler) -> Self {
        Self {
            matrix_backend: ReactorBvpMatrixBackend::Sparse,
            execution_backend: ReactorBvpExecutionBackend::Aot(ReactorBvpAotConfig::for_compiler(
                compiler,
            )),
            ..Self::default_lambdify()
        }
    }

    /// Return a config with another matrix backend.
    pub fn with_matrix_backend(mut self, matrix_backend: ReactorBvpMatrixBackend) -> Self {
        self.matrix_backend = matrix_backend;
        self
    }

    /// Return a config with another symbolic assembly backend.
    pub fn with_symbolic_backend(mut self, symbolic_backend: ReactorBvpSymbolicBackend) -> Self {
        self.symbolic_backend = symbolic_backend;
        self
    }

    /// Return a config with another generated execution backend.
    pub fn with_execution_backend(mut self, execution_backend: ReactorBvpExecutionBackend) -> Self {
        self.execution_backend = execution_backend;
        self
    }

    /// Return a config with updated AOT settings, promoting lambdify to banded AOT when needed.
    pub fn with_aot_config(mut self, aot_config: ReactorBvpAotConfig) -> Self {
        self.execution_backend = ReactorBvpExecutionBackend::Aot(aot_config);
        self
    }

    /// Mutate the current AOT config, creating the default AOT branch when absent.
    pub fn map_aot_config(
        mut self,
        f: impl FnOnce(ReactorBvpAotConfig) -> ReactorBvpAotConfig,
    ) -> Self {
        let current = match self.execution_backend {
            ReactorBvpExecutionBackend::Aot(config) => config,
            ReactorBvpExecutionBackend::Lambdify => ReactorBvpAotConfig::default(),
        };
        self.execution_backend = ReactorBvpExecutionBackend::Aot(f(current));
        self
    }

    /// Parse high-level generated backend aliases used by task documents.
    pub fn from_generated_backend_name(value: &str) -> Result<Self, String> {
        let config = Self::default_lambdify();
        match normalize_key(value).as_str() {
            "banded_lambdify" | "lambdify_banded" | "lambdify" => Ok(config),
            "sparse_lambdify" | "lambdify_sparse" => Ok(Self::sparse_lambdify()),
            "banded_aot" | "aot_banded" => Ok(Self::banded_aot(ReactorBvpAotCompiler::CTcc)),
            "sparse_aot" | "aot_sparse" => Ok(Self::sparse_aot(ReactorBvpAotCompiler::CTcc)),
            "banded_aot_tcc" | "aot_banded_tcc" => {
                Ok(Self::banded_aot(ReactorBvpAotCompiler::CTcc))
            }
            "banded_aot_gcc" | "aot_banded_gcc" => {
                Ok(Self::banded_aot(ReactorBvpAotCompiler::CGcc))
            }
            "banded_aot_zig" | "aot_banded_zig" => Ok(Self::banded_aot(ReactorBvpAotCompiler::Zig)),
            "sparse_aot_tcc" | "aot_sparse_tcc" => {
                Ok(Self::sparse_aot(ReactorBvpAotCompiler::CTcc))
            }
            "sparse_aot_gcc" | "aot_sparse_gcc" => {
                Ok(Self::sparse_aot(ReactorBvpAotCompiler::CGcc))
            }
            "sparse_aot_zig" | "aot_sparse_zig" => Ok(Self::sparse_aot(ReactorBvpAotCompiler::Zig)),
            other => Err(format!(
                "Unknown reactor BVP generated backend `{other}`. Expected `banded_lambdify`, `sparse_lambdify`, `banded_aot_*`, or `sparse_aot_*`"
            )),
        }
    }

    /// Convert the KiThe facade into RustedSciThe damped-solver options.
    pub fn to_rusted_options(self) -> Result<DampedSolverOptions, String> {
        let options = match self.execution_backend {
            ReactorBvpExecutionBackend::Lambdify => self.to_lambdify_options(),
            ReactorBvpExecutionBackend::Aot(config) => self.to_aot_options(config)?,
        };
        Ok(options)
    }

    fn to_lambdify_options(self) -> DampedSolverOptions {
        match self.matrix_backend {
            ReactorBvpMatrixBackend::Banded => DampedSolverOptions::banded_damped()
                .with_banded_lambdify()
                .with_symbolic_assembly_backend(self.symbolic_backend.to_rusted()),
            ReactorBvpMatrixBackend::Sparse => {
                let generated = GeneratedBackendConfig::sparse_defaults()
                    .with_backend_policy_override(Some(BackendSelectionPolicy::LambdifyOnly))
                    .with_symbolic_assembly_backend(self.symbolic_backend.to_rusted());
                DampedSolverOptions::sparse_damped().with_generated_backend_config(generated)
            }
        }
    }

    fn to_aot_options(
        self,
        aot_config: ReactorBvpAotConfig,
    ) -> Result<DampedSolverOptions, String> {
        let generated = base_aot_generated_config(self.matrix_backend, aot_config.compiler)
            .with_symbolic_assembly_backend(self.symbolic_backend.to_rusted())
            .with_aot_build_policy(aot_config.build_policy.to_rusted(aot_config.build_profile))
            .with_aot_execution_policy(aot_config.execution_policy.to_rusted())
            .with_aot_chunking_policy(aot_chunking_policy(aot_config));

        let generated = match aot_config.compile_preset {
            ReactorBvpAotCompilePreset::Production => generated.with_aot_compile_production(),
            ReactorBvpAotCompilePreset::FastBuild => generated.with_aot_compile_fast_build(),
            ReactorBvpAotCompilePreset::DevFastest => generated.with_aot_compile_dev_fastest(),
        };

        let options = match self.matrix_backend {
            ReactorBvpMatrixBackend::Banded => DampedSolverOptions::banded_damped(),
            ReactorBvpMatrixBackend::Sparse => DampedSolverOptions::sparse_damped(),
        }
        .with_generated_backend_config(generated);
        Ok(options)
    }
}

fn base_aot_generated_config(
    matrix_backend: ReactorBvpMatrixBackend,
    compiler: ReactorBvpAotCompiler,
) -> GeneratedBackendConfig {
    match (matrix_backend, compiler) {
        (ReactorBvpMatrixBackend::Banded, ReactorBvpAotCompiler::CTcc) => {
            GeneratedBackendConfig::banded_atomview_build_if_missing_release_tcc()
        }
        (ReactorBvpMatrixBackend::Banded, ReactorBvpAotCompiler::CGcc) => {
            GeneratedBackendConfig::banded_atomview_build_if_missing_release_gcc()
        }
        (ReactorBvpMatrixBackend::Banded, ReactorBvpAotCompiler::Zig) => {
            GeneratedBackendConfig::banded_atomview_build_if_missing_release_zig()
        }
        (ReactorBvpMatrixBackend::Sparse, ReactorBvpAotCompiler::CTcc) => {
            GeneratedBackendConfig::sparse_atomview_build_if_missing_release_tcc()
        }
        (ReactorBvpMatrixBackend::Sparse, ReactorBvpAotCompiler::CGcc) => {
            GeneratedBackendConfig::sparse_atomview_build_if_missing_release_gcc()
        }
        (ReactorBvpMatrixBackend::Sparse, ReactorBvpAotCompiler::Zig) => {
            GeneratedBackendConfig::sparse_atomview_build_if_missing_release_zig()
        }
    }
}

fn aot_chunking_policy(config: ReactorBvpAotConfig) -> AotChunkingPolicy {
    let residual = config
        .residual_target_chunks
        .map(|target_chunks| ResidualChunkingStrategy::ByTargetChunkCount { target_chunks });
    let jacobian = config
        .jacobian_target_chunks
        .map(|target_chunks| SparseChunkingStrategy::ByTargetChunkCount { target_chunks });
    AotChunkingPolicy::with_parts(residual, jacobian)
}

fn normalize_key(value: &str) -> String {
    value.trim().to_ascii_lowercase().replace('-', "_")
}
