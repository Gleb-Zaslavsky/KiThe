//! Typed editor model for the chemical-equilibrium GUI.
//!
//! This module deliberately stops before the repository and solver boundary.
//! The document stores editable values and policies; validation produces a
//! separate SI-valued model that a later request builder can pass to
//! `Thermodynamics::ChemEquilibrium::prelude`.
//!
//! Keeping this boundary pure gives the GUI three useful properties:
//! no database access while typing, deterministic invalidation of derived
//! candidate previews, and no possibility of accidentally solving with a
//! half-edited numeric field.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashSet};
use std::fmt;

/// Current on-disk schema version for an equilibrium GUI document.
pub const EQUILIBRIUM_GUI_SCHEMA_VERSION: u32 = 1;

/// A complete editable equilibrium document.
///
/// Runtime workers, repository handles, candidate previews, plots, and result
/// snapshots are intentionally not fields of this value. Saving and loading a
/// document therefore cannot resurrect a stale solve.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumGuiDocument {
    /// Schema version used for future migrations.
    pub schema_version: u32,
    /// Editable configuration owned by the document.
    pub config: EquilibriumGuiConfig,
}

impl Default for EquilibriumGuiDocument {
    fn default() -> Self {
        Self {
            schema_version: EQUILIBRIUM_GUI_SCHEMA_VERSION,
            config: EquilibriumGuiConfig::default(),
        }
    }
}

impl EquilibriumGuiDocument {
    /// Creates the current schema version with default, runnable `P,T` input.
    pub fn new() -> Self {
        Self::default()
    }

    /// Creates the simple explicit ideal-gas preset.
    ///
    /// This is intentionally only a convenience constructor. The returned
    /// document uses the same `ExplicitPhases` representation and therefore
    /// follows exactly the same validation and request-building path as the
    /// advanced phase editor.
    pub fn simple_ideal_gas<I, S>(substances: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let mut components = substances
            .into_iter()
            .map(|substance| ComponentDraft {
                substance: substance.into(),
                initial_moles: "1.0".into(),
                source_library: None,
            })
            .collect::<Vec<_>>();
        if components.is_empty() {
            components.push(ComponentDraft::default());
        }
        Self {
            schema_version: EQUILIBRIUM_GUI_SCHEMA_VERSION,
            config: EquilibriumGuiConfig {
                inventory: EquilibriumInventoryDraft::ExplicitPhases {
                    phases: vec![PhaseDraft {
                        id: "gas".into(),
                        physical_state: GuiPhysicalState::Gas,
                        model: GuiPhaseModel::IdealGas,
                        components,
                    }],
                },
                ..EquilibriumGuiConfig::default()
            },
        }
    }

    /// Rejects documents from a future schema before they enter the editor.
    pub fn validate_schema(&self) -> Result<(), EquilibriumGuiValidationReport> {
        if self.schema_version != EQUILIBRIUM_GUI_SCHEMA_VERSION {
            return Err(EquilibriumGuiValidationReport::single(
                "schema_version",
                ValidationIssueKind::UnsupportedSchema,
                format!(
                    "document schema {} is not supported; expected {}",
                    self.schema_version, EQUILIBRIUM_GUI_SCHEMA_VERSION
                ),
            ));
        }
        Ok(())
    }

    /// Validates fields needed for editing and candidate preview.
    pub fn validate(
        &self,
    ) -> Result<ValidatedEquilibriumGuiConfig, EquilibriumGuiValidationReport> {
        self.validate_schema()?;
        self.config.validate(false)
    }

    /// Validates that the document can be sent to a production solve.
    ///
    /// Element mode is allowed to exist without a preview while the user is
    /// editing. A solve requires that the preview has been confirmed into
    /// explicit phase/component assignments.
    pub fn validate_for_run(
        &self,
    ) -> Result<ValidatedEquilibriumGuiConfig, EquilibriumGuiValidationReport> {
        self.validate_schema()?;
        self.config.validate(true)
    }

    /// Stable JSON representation used by document lifecycle code and tests.
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Loads and schema-checks a saved document.
    pub fn from_json(json: &str) -> Result<Self, EquilibriumGuiDocumentError> {
        let document: Self =
            serde_json::from_str(json).map_err(EquilibriumGuiDocumentError::Deserialize)?;
        document
            .validate_schema()
            .map_err(EquilibriumGuiDocumentError::InvalidSchema)?;
        Ok(document)
    }
}

/// Editable configuration grouped by domain policy.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumGuiConfig {
    pub problem: EquilibriumProblemDraft,
    pub inventory: EquilibriumInventoryDraft,
    pub lookup: EquilibriumLookupDraft,
    pub phase_mode: EquilibriumPhaseModeDraft,
    pub solver: EquilibriumSolverConfigDraft,
    pub diagnostics: EquilibriumDiagnosticsDraft,
    pub postprocessing: EquilibriumPostprocessingDraft,
}

impl Default for EquilibriumGuiConfig {
    fn default() -> Self {
        Self {
            problem: EquilibriumProblemDraft::default(),
            inventory: EquilibriumInventoryDraft::default(),
            lookup: EquilibriumLookupDraft::default(),
            phase_mode: EquilibriumPhaseModeDraft::default(),
            solver: EquilibriumSolverConfigDraft::default(),
            diagnostics: EquilibriumDiagnosticsDraft::default(),
            postprocessing: EquilibriumPostprocessingDraft::default(),
        }
    }
}

impl EquilibriumGuiConfig {
    fn validate(
        &self,
        for_run: bool,
    ) -> Result<ValidatedEquilibriumGuiConfig, EquilibriumGuiValidationReport> {
        let mut report = EquilibriumGuiValidationReport::default();
        let problem = self.problem.validate(&mut report);
        let inventory = self.inventory.validate(&mut report, for_run);
        let lookup = self.lookup.validate(&mut report);
        let phase_mode = self.phase_mode.validate(&mut report);
        let solver = self.solver.validate(&mut report);
        let diagnostics = self.diagnostics.validate(&mut report);
        let postprocessing = self.postprocessing.validate(&mut report);

        if report.has_errors() {
            return Err(report);
        }

        Ok(ValidatedEquilibriumGuiConfig {
            problem: problem.expect("validated problem is present when report has no errors"),
            inventory: inventory.expect("validated inventory is present when report has no errors"),
            lookup: lookup.expect("validated lookup is present when report has no errors"),
            phase_mode: phase_mode
                .expect("validated phase mode is present when report has no errors"),
            solver: solver.expect("validated solver is present when report has no errors"),
            diagnostics: diagnostics
                .expect("validated diagnostics are present when report has no errors"),
            postprocessing: postprocessing
                .expect("validated postprocessing is present when report has no errors"),
        })
    }
}

/// Fixed-pressure thermodynamic problem draft.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum EquilibriumProblemDraft {
    /// Production-supported ideal-system constraint.
    FixedPt {
        pressure_pa: String,
        reference_pressure_pa: String,
        temperature: TemperatureDraft,
    },
    /// Fixed-pressure, fixed-total-enthalpy problem.
    ///
    /// `target_enthalpy_j` is extensive total enthalpy in joules. The outer
    /// temperature solve is constrained to the common thermochemistry domain
    /// described by `temperature_bounds`; `seed_temperature_k` is only the
    /// initial scalar-solver estimate.
    FixedPh {
        pressure_pa: String,
        reference_pressure_pa: String,
        #[serde(alias = "enthalpy_j_per_mol")]
        target_enthalpy_j: String,
        #[serde(default)]
        temperature_bounds: PhTemperatureBoundsDraft,
    },
}

impl Default for EquilibriumProblemDraft {
    fn default() -> Self {
        Self::FixedPt {
            pressure_pa: "101325".into(),
            reference_pressure_pa: "101325".into(),
            temperature: TemperatureDraft::default(),
        }
    }
}

/// Point or ordered temperature grid draft.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum TemperatureDraft {
    Point {
        temperature_k: String,
    },
    Range {
        start_k: String,
        end_k: String,
        point_count: String,
    },
}

impl Default for TemperatureDraft {
    fn default() -> Self {
        Self::Point {
            temperature_k: "1000".into(),
        }
    }
}

/// Editable temperature bracket and initial estimate for a `P,H` solve.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PhTemperatureBoundsDraft {
    pub lower_k: String,
    pub upper_k: String,
    pub seed_k: String,
}

impl Default for PhTemperatureBoundsDraft {
    fn default() -> Self {
        Self {
            lower_k: "300".into(),
            upper_k: "3000".into(),
            seed_k: "1000".into(),
        }
    }
}

/// Inventory input mode. Element mode becomes runnable only after candidate
/// selection has produced explicit phase assignments.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum EquilibriumInventoryDraft {
    ExplicitPhases {
        phases: Vec<PhaseDraft>,
    },
    ElementCandidates {
        elements: Vec<String>,
        candidate_policy: CandidatePolicyDraft,
        assignments: Vec<PhaseDraft>,
    },
}

impl Default for EquilibriumInventoryDraft {
    fn default() -> Self {
        Self::ExplicitPhases {
            phases: vec![PhaseDraft::default()],
        }
    }
}

/// One editable phase and its initial phase-qualified components.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PhaseDraft {
    /// Empty ID means the anonymous single-phase identity.
    pub id: String,
    pub physical_state: GuiPhysicalState,
    pub model: GuiPhaseModel,
    pub components: Vec<ComponentDraft>,
}

impl Default for PhaseDraft {
    fn default() -> Self {
        Self {
            id: "gas".into(),
            physical_state: GuiPhysicalState::Gas,
            model: GuiPhaseModel::IdealGas,
            components: vec![ComponentDraft::default()],
        }
    }
}

/// One exact catalog record and its initial physical amount.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ComponentDraft {
    pub substance: String,
    pub initial_moles: String,
    /// Optional catalog provenance pinned by the element-candidate preview.
    ///
    /// Explicitly typed components may leave this unset and use the normal
    /// lookup policy. Candidate-driven components carry the library selected
    /// by the immutable repository report so a later resolve cannot silently
    /// substitute a record from another library.
    #[serde(default)]
    pub source_library: Option<String>,
}

impl Default for ComponentDraft {
    fn default() -> Self {
        Self {
            substance: "H2O".into(),
            initial_moles: "1.0".into(),
            source_library: None,
        }
    }
}

/// Catalog filtering policy used by element candidate preview.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CandidatePolicyDraft {
    pub element_mode: GuiElementSearchMode,
    pub physical_states: Vec<GuiPhysicalState>,
    pub temperature_lower_k: String,
    pub temperature_upper_k: String,
    pub max_candidates: String,
}

impl Default for CandidatePolicyDraft {
    fn default() -> Self {
        Self {
            element_mode: GuiElementSearchMode::SubsetOf,
            physical_states: vec![GuiPhysicalState::Gas],
            temperature_lower_k: "300".into(),
            temperature_upper_k: "3000".into(),
            max_candidates: "100".into(),
        }
    }
}

/// Lookup policy shared by every phase in one document.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum EquilibriumLookupDraft {
    Default,
    Explicit {
        priority_libraries: Vec<String>,
        permitted_libraries: Vec<String>,
        explicit_search_instructions: BTreeMap<String, String>,
        search_in_nist: bool,
    },
}

impl Default for EquilibriumLookupDraft {
    fn default() -> Self {
        Self::Default
    }
}

/// Fixed declared phases or bounded active-set phase control.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum EquilibriumPhaseModeDraft {
    FixedDeclared,
    Bounded {
        phase_epsilon: String,
        dg_create: String,
        dg_keep: String,
        max_phase_iterations: String,
    },
}

impl Default for EquilibriumPhaseModeDraft {
    fn default() -> Self {
        Self::FixedDeclared
    }
}

/// Concrete nonlinear backend selected by the GUI.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GuiSolverBackend {
    RstLm,
    RstMinpackLm,
    RstNielsenLm,
    RstTrustRegionLm,
    RstPowellDogleg,
    RstDampedNewton,
    LegacyLm,
    LegacyNr,
    LegacyTr,
}

impl GuiSolverBackend {
    /// Stable order used by the GUI selector, serialization stories, and
    /// backend-matrix tests. Keep this list exhaustive when the engine adds a
    /// supported fallback.
    pub const ALL: [Self; 9] = [
        Self::RstLm,
        Self::RstMinpackLm,
        Self::RstNielsenLm,
        Self::RstTrustRegionLm,
        Self::RstPowellDogleg,
        Self::RstDampedNewton,
        Self::LegacyLm,
        Self::LegacyNr,
        Self::LegacyTr,
    ];

    pub const fn label(self) -> &'static str {
        match self {
            Self::RstLm => "RST Levenberg-Marquardt",
            Self::RstMinpackLm => "RST Minpack Levenberg-Marquardt",
            Self::RstNielsenLm => "RST Nielsen Levenberg-Marquardt",
            Self::RstTrustRegionLm => "RST Trust-region LM",
            Self::RstPowellDogleg => "RST Powell dogleg",
            Self::RstDampedNewton => "RST damped Newton",
            Self::LegacyLm => "Legacy LM fallback",
            Self::LegacyNr => "Legacy Newton-Raphson fallback",
            Self::LegacyTr => "Legacy trust-region fallback",
        }
    }
}

/// Production default, one explicit backend, or a user-owned ordered cascade.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum EquilibriumSolverDraft {
    ProductionDefault,
    SingleBackend { backend: GuiSolverBackend },
    CustomCascade { backends: Vec<GuiSolverBackend> },
}

impl Default for EquilibriumSolverDraft {
    fn default() -> Self {
        Self::ProductionDefault
    }
}

/// Optional advanced solver controls. Empty strings mean engine defaults.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumSolverOverridesDraft {
    pub tolerance: String,
    pub max_iterations: String,
    pub scaling_enabled: bool,
    #[serde(default)]
    pub trace_seed_policy: Option<GuiTraceSeedPolicyDraft>,
    /// Optional explicit limits for a custom or production backend cascade.
    #[serde(default)]
    pub cascade_budget: Option<GuiSolverCascadeBudgetDraft>,
}

impl Default for EquilibriumSolverOverridesDraft {
    fn default() -> Self {
        Self {
            tolerance: String::new(),
            max_iterations: String::new(),
            // Match the equilibrium engine's stable default. Scaling remains
            // available as an explicit advanced option because bounded
            // multiphase solves can be more sensitive to that transformation.
            scaling_enabled: false,
            trace_seed_policy: None,
            cascade_budget: None,
        }
    }
}

/// Serializable editor form for the engine's deterministic cascade budget.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GuiSolverCascadeBudgetDraft {
    pub max_attempts: String,
    pub max_iterations_per_attempt: String,
    pub max_total_iterations: String,
}

/// Validated cascade limits handed to the request builder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ValidatedSolverCascadeBudget {
    pub max_attempts: usize,
    pub max_iterations_per_attempt: usize,
    pub max_total_iterations: usize,
}

/// Optional policy for positive log-mole coordinates. `None` keeps the
/// canonical engine default and therefore does not duplicate solver policy in
/// the GUI document.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum GuiTraceSeedPolicyDraft {
    Absolute {
        floor: String,
    },
    RelativeToLargestInitialMole {
        fraction: String,
        minimum_floor: String,
    },
}

/// Solver policy plus optional validated overrides.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumSolverConfigDraft {
    pub selection: EquilibriumSolverDraft,
    pub overrides: EquilibriumSolverOverridesDraft,
}

impl Default for EquilibriumSolverConfigDraft {
    fn default() -> Self {
        Self {
            selection: EquilibriumSolverDraft::default(),
            overrides: EquilibriumSolverOverridesDraft::default(),
        }
    }
}

/// A compact diagnostics policy; these flags control retained report detail,
/// not process-global logger configuration.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumDiagnosticsDraft {
    pub collect_timing: bool,
    pub retain_backend_attempts: bool,
    pub retain_conservation_report: bool,
    pub retain_phase_transitions: bool,
    pub keq_validation: GuiKeqValidationMode,
}

impl Default for EquilibriumDiagnosticsDraft {
    fn default() -> Self {
        Self {
            collect_timing: false,
            retain_backend_attempts: true,
            retain_conservation_report: true,
            retain_phase_transitions: true,
            keq_validation: GuiKeqValidationMode::WhenApplicable,
        }
    }
}

/// Independent K_eq validation policy exposed without importing engine
/// internals. The request builder will map this to the production enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiKeqValidationMode {
    Off,
    WhenApplicable,
    Required,
}

/// Display and resampling policy for accepted range results.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumPostprocessingDraft {
    pub plot_target: GuiPlotTarget,
    pub result_basis: GuiResultBasis,
    pub resampling: GuiResamplingDraft,
    #[serde(default)]
    pub y_scale: GuiPlotScale,
}

impl Default for EquilibriumPostprocessingDraft {
    fn default() -> Self {
        Self {
            plot_target: GuiPlotTarget::None,
            result_basis: GuiResultBasis::ComponentMoles,
            resampling: GuiResamplingDraft::None,
            y_scale: GuiPlotScale::Linear,
        }
    }
}

/// Which existing plot surface should receive a result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiPlotTarget {
    None,
    Embedded,
    KiThePlot,
    Both,
}

/// Quantity represented by a result plot/table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiResultBasis {
    ComponentMoles,
    MoleFractions,
    PhaseTotals,
}

/// Display-only ordinate scale shared by both equilibrium plot adapters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiPlotScale {
    Linear,
    Log10,
}

impl Default for GuiPlotScale {
    fn default() -> Self {
        Self::Linear
    }
}

/// Optional display-only range resampling.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum GuiResamplingDraft {
    None,
    Pchip {
        output_points: String,
        #[serde(default)]
        interpolation_space: GuiInterpolationSpace,
        #[serde(default)]
        clamp: bool,
    },
}

impl Default for GuiResamplingDraft {
    fn default() -> Self {
        Self::None
    }
}

/// Coordinate space for display-only PCHIP interpolation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiInterpolationSpace {
    Linear,
    Log,
}

impl Default for GuiInterpolationSpace {
    fn default() -> Self {
        Self::Linear
    }
}

/// Physical state chosen in the GUI. This is converted to the engine state at
/// the request boundary instead of deriving it from a record name.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GuiPhysicalState {
    Gas,
    Liquid,
    Solid,
    Condensed,
}

/// Activity model currently supported by the phase engine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiPhaseModel {
    IdealGas,
    PureCondensed,
}

/// Element matching policy used by candidate selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GuiElementSearchMode {
    Exact,
    SubsetOf,
}

/// Validated SI-valued model handed to the request builder.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedEquilibriumGuiConfig {
    pub problem: ValidatedProblem,
    pub inventory: ValidatedInventory,
    pub lookup: ValidatedLookup,
    pub phase_mode: ValidatedPhaseMode,
    pub solver: ValidatedSolver,
    pub diagnostics: ValidatedDiagnostics,
    pub postprocessing: ValidatedPostprocessing,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedProblem {
    FixedPt {
        pressure_pa: f64,
        reference_pressure_pa: f64,
        temperature: ValidatedTemperature,
    },
    FixedPh {
        pressure_pa: f64,
        reference_pressure_pa: f64,
        target_enthalpy_j: f64,
        temperature_bounds: ValidatedPhTemperatureBounds,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ValidatedPhTemperatureBounds {
    pub lower_k: f64,
    pub upper_k: f64,
    pub seed_k: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedTemperature {
    Point(f64),
    Range {
        start_k: f64,
        end_k: f64,
        point_count: usize,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedInventory {
    ExplicitPhases(Vec<ValidatedPhase>),
    ElementCandidates {
        elements: Vec<String>,
        candidate_policy: ValidatedCandidatePolicy,
        assignments: Vec<ValidatedPhase>,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedPhase {
    pub id: Option<String>,
    pub physical_state: GuiPhysicalState,
    pub model: GuiPhaseModel,
    pub components: Vec<ValidatedComponent>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedComponent {
    pub substance: String,
    pub initial_moles: f64,
    pub source_library: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedCandidatePolicy {
    pub element_mode: GuiElementSearchMode,
    pub physical_states: Vec<GuiPhysicalState>,
    pub temperature_lower_k: f64,
    pub temperature_upper_k: f64,
    pub max_candidates: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedLookup {
    Default,
    Explicit {
        priority_libraries: Vec<String>,
        permitted_libraries: Vec<String>,
        explicit_search_instructions: BTreeMap<String, String>,
        search_in_nist: bool,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedPhaseMode {
    FixedDeclared,
    Bounded {
        phase_epsilon: f64,
        dg_create: f64,
        dg_keep: f64,
        max_phase_iterations: usize,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedSolver {
    pub selection: EquilibriumSolverDraft,
    pub tolerance: Option<f64>,
    pub max_iterations: Option<usize>,
    pub scaling_enabled: bool,
    pub trace_seed_policy: Option<ValidatedTraceSeedPolicy>,
    pub cascade_budget: Option<ValidatedSolverCascadeBudget>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ValidatedTraceSeedPolicy {
    Absolute { floor: f64 },
    RelativeToLargestInitialMole { fraction: f64, minimum_floor: f64 },
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedDiagnostics {
    pub collect_timing: bool,
    pub retain_backend_attempts: bool,
    pub retain_conservation_report: bool,
    pub retain_phase_transitions: bool,
    pub keq_validation: GuiKeqValidationMode,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedPostprocessing {
    pub plot_target: GuiPlotTarget,
    pub result_basis: GuiResultBasis,
    pub resampling: Option<usize>,
    pub y_scale: GuiPlotScale,
}

/// Severity/category of one pure validation finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ValidationIssueKind {
    MissingValue,
    InvalidValue,
    Duplicate,
    Incompatible,
    FeatureUnavailable,
    UnsupportedSchema,
}

/// Field-addressable validation finding suitable for inline GUI display.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ValidationIssue {
    pub field: String,
    pub kind: ValidationIssueKind,
    pub message: String,
}

/// Aggregate pure validation result.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EquilibriumGuiValidationReport {
    pub issues: Vec<ValidationIssue>,
}

impl EquilibriumGuiValidationReport {
    fn single(
        field: impl Into<String>,
        kind: ValidationIssueKind,
        message: impl Into<String>,
    ) -> Self {
        let mut report = Self::default();
        report.push(field, kind, message);
        report
    }

    fn push(
        &mut self,
        field: impl Into<String>,
        kind: ValidationIssueKind,
        message: impl Into<String>,
    ) {
        self.issues.push(ValidationIssue {
            field: field.into(),
            kind,
            message: message.into(),
        });
    }

    pub fn has_errors(&self) -> bool {
        !self.issues.is_empty()
    }

    pub fn contains_kind(&self, kind: ValidationIssueKind) -> bool {
        self.issues.iter().any(|issue| issue.kind == kind)
    }
}

/// Errors raised while loading a GUI document.
#[derive(Debug)]
pub enum EquilibriumGuiDocumentError {
    Deserialize(serde_json::Error),
    InvalidSchema(EquilibriumGuiValidationReport),
}

impl fmt::Display for EquilibriumGuiDocumentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Deserialize(error) => {
                write!(f, "failed to deserialize equilibrium document: {error}")
            }
            Self::InvalidSchema(report) => {
                write!(
                    f,
                    "invalid equilibrium document schema: {:?}",
                    report.issues
                )
            }
        }
    }
}

impl std::error::Error for EquilibriumGuiDocumentError {}

impl EquilibriumProblemDraft {
    fn validate(&self, report: &mut EquilibriumGuiValidationReport) -> Option<ValidatedProblem> {
        match self {
            Self::FixedPt {
                pressure_pa,
                reference_pressure_pa,
                temperature,
            } => {
                let pressure = parse_positive(pressure_pa, "problem.pressure_pa", report);
                let reference = parse_positive(
                    reference_pressure_pa,
                    "problem.reference_pressure_pa",
                    report,
                );
                let temperature = temperature.validate(report);
                match (pressure, reference, temperature) {
                    (Some(pressure_pa), Some(reference_pressure_pa), Some(temperature)) => {
                        Some(ValidatedProblem::FixedPt {
                            pressure_pa,
                            reference_pressure_pa,
                            temperature,
                        })
                    }
                    _ => None,
                }
            }
            Self::FixedPh {
                pressure_pa,
                reference_pressure_pa,
                target_enthalpy_j,
                temperature_bounds,
            } => {
                let pressure = parse_positive(pressure_pa, "problem.pressure_pa", report);
                let reference = parse_positive(
                    reference_pressure_pa,
                    "problem.reference_pressure_pa",
                    report,
                );
                let enthalpy = parse_finite(target_enthalpy_j, "problem.target_enthalpy_j", report);
                let lower = parse_positive(
                    &temperature_bounds.lower_k,
                    "problem.temperature_bounds.lower_k",
                    report,
                );
                let upper = parse_positive(
                    &temperature_bounds.upper_k,
                    "problem.temperature_bounds.upper_k",
                    report,
                );
                let seed = parse_positive(
                    &temperature_bounds.seed_k,
                    "problem.temperature_bounds.seed_k",
                    report,
                );
                let ordered = matches!((lower, upper), (Some(lower), Some(upper)) if lower < upper);
                if matches!((lower, upper), (Some(lower), Some(upper)) if lower >= upper) {
                    report.push(
                        "problem.temperature_bounds",
                        ValidationIssueKind::InvalidValue,
                        "P,H lower temperature bound must be strictly below upper bound",
                    );
                }
                if let (Some(lower), Some(upper), Some(seed)) = (lower, upper, seed) {
                    if !(lower..=upper).contains(&seed) {
                        report.push(
                            "problem.temperature_bounds.seed_k",
                            ValidationIssueKind::InvalidValue,
                            "P,H seed temperature must lie inside the temperature bounds",
                        );
                    }
                }
                match (pressure, reference, enthalpy, lower, upper, seed, ordered) {
                    (
                        Some(pressure_pa),
                        Some(reference_pressure_pa),
                        Some(target_enthalpy_j),
                        Some(lower_k),
                        Some(upper_k),
                        Some(seed_k),
                        true,
                    ) => Some(ValidatedProblem::FixedPh {
                        pressure_pa,
                        reference_pressure_pa,
                        target_enthalpy_j,
                        temperature_bounds: ValidatedPhTemperatureBounds {
                            lower_k,
                            upper_k,
                            seed_k,
                        },
                    }),
                    _ => None,
                }
            }
        }
    }
}

impl TemperatureDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedTemperature> {
        match self {
            Self::Point { temperature_k } => parse_positive(temperature_k, "temperature_k", report)
                .map(ValidatedTemperature::Point),
            Self::Range {
                start_k,
                end_k,
                point_count,
            } => {
                let start = parse_positive(start_k, "temperature.start_k", report);
                let end = parse_positive(end_k, "temperature.end_k", report);
                let count = parse_usize_at_least(point_count, "temperature.point_count", 2, report);
                let equal_endpoints =
                    matches!((start, end), (Some(start), Some(end)) if start == end);
                if equal_endpoints {
                    report.push(
                        "temperature",
                        ValidationIssueKind::InvalidValue,
                        "temperature range endpoints must differ",
                    );
                }
                match (start, end, count, equal_endpoints) {
                    (Some(start_k), Some(end_k), Some(point_count), false) => {
                        Some(ValidatedTemperature::Range {
                            start_k,
                            end_k,
                            point_count,
                        })
                    }
                    _ => None,
                }
            }
        }
    }
}

impl EquilibriumInventoryDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
        for_run: bool,
    ) -> Option<ValidatedInventory> {
        match self {
            Self::ExplicitPhases { phases } => {
                validate_phases(phases, report).map(ValidatedInventory::ExplicitPhases)
            }
            Self::ElementCandidates {
                elements,
                candidate_policy,
                assignments,
            } => {
                let mut unique_elements = Vec::with_capacity(elements.len());
                let mut seen = HashSet::new();
                for (index, element) in elements.iter().enumerate() {
                    let element = element.trim();
                    if element.is_empty() {
                        report.push(
                            format!("inventory.elements[{index}]"),
                            ValidationIssueKind::MissingValue,
                            "element symbol must not be empty",
                        );
                    } else if !seen.insert(element.to_string()) {
                        report.push(
                            format!("inventory.elements[{index}]"),
                            ValidationIssueKind::Duplicate,
                            format!("element '{element}' is listed more than once"),
                        );
                    } else {
                        unique_elements.push(element.to_string());
                    }
                }
                if unique_elements.is_empty() {
                    report.push(
                        "inventory.elements",
                        ValidationIssueKind::MissingValue,
                        "at least one element is required",
                    );
                }
                let policy = candidate_policy.validate(report);
                let validated_assignments = if assignments.is_empty() {
                    Vec::new()
                } else {
                    validate_phases(assignments, report).unwrap_or_default()
                };
                if for_run && assignments.is_empty() {
                    report.push(
                        "inventory.assignments",
                        ValidationIssueKind::MissingValue,
                        "confirm candidate selection and assign at least one phase before running",
                    );
                }
                policy.map(|candidate_policy| ValidatedInventory::ElementCandidates {
                    elements: unique_elements,
                    candidate_policy,
                    assignments: validated_assignments,
                })
            }
        }
    }
}

fn validate_phases(
    phases: &[PhaseDraft],
    report: &mut EquilibriumGuiValidationReport,
) -> Option<Vec<ValidatedPhase>> {
    if phases.is_empty() {
        report.push(
            "inventory.phases",
            ValidationIssueKind::MissingValue,
            "at least one phase is required",
        );
        return None;
    }

    let mut seen_phase_ids = HashSet::new();
    let mut validated = Vec::with_capacity(phases.len());
    let mut any_positive_moles = false;
    let mut seen_components = HashSet::new();
    for (phase_index, phase) in phases.iter().enumerate() {
        let id = phase.id.trim();
        let phase_key = if id.is_empty() { "<anonymous>" } else { id };
        if !seen_phase_ids.insert(phase_key.to_string()) {
            report.push(
                format!("inventory.phases[{phase_index}].id"),
                ValidationIssueKind::Duplicate,
                format!("phase '{phase_key}' is declared more than once"),
            );
        }
        if !phase_model_matches_state(phase.physical_state, phase.model) {
            report.push(
                format!("inventory.phases[{phase_index}].model"),
                ValidationIssueKind::Incompatible,
                "the selected activity model is incompatible with the physical state",
            );
        }
        if phase.components.is_empty() {
            report.push(
                format!("inventory.phases[{phase_index}].components"),
                ValidationIssueKind::MissingValue,
                "each phase needs at least one component",
            );
        }
        let mut components = Vec::with_capacity(phase.components.len());
        for (component_index, component) in phase.components.iter().enumerate() {
            let substance = component.substance.trim();
            if substance.is_empty() {
                report.push(
                    format!(
                        "inventory.phases[{phase_index}].components[{component_index}].substance"
                    ),
                    ValidationIssueKind::MissingValue,
                    "substance name must not be empty",
                );
                continue;
            }
            let qualified = format!("{phase_key}::{substance}");
            if !seen_components.insert(qualified.clone()) {
                report.push(
                    format!("inventory.phases[{phase_index}].components[{component_index}]"),
                    ValidationIssueKind::Duplicate,
                    format!("phase-qualified component '{qualified}' is duplicated"),
                );
            }
            if let Some(moles) = parse_non_negative(
                &component.initial_moles,
                &format!(
                    "inventory.phases[{phase_index}].components[{component_index}].initial_moles"
                ),
                report,
            ) {
                any_positive_moles |= moles > 0.0;
                components.push(ValidatedComponent {
                    substance: substance.to_string(),
                    initial_moles: moles,
                    source_library: validate_source_library(
                        component.source_library.as_deref(),
                        phase_index,
                        component_index,
                        report,
                    ),
                });
            }
        }
        validated.push(ValidatedPhase {
            id: (!id.is_empty()).then(|| id.to_string()),
            physical_state: phase.physical_state,
            model: phase.model,
            components,
        });
    }
    if !any_positive_moles {
        report.push(
            "inventory",
            ValidationIssueKind::InvalidValue,
            "at least one initial component amount must be positive",
        );
    }
    Some(validated)
}

fn validate_source_library(
    source_library: Option<&str>,
    phase_index: usize,
    component_index: usize,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<String> {
    let Some(source_library) = source_library else {
        return None;
    };
    let source_library = source_library.trim();
    if source_library.is_empty() {
        report.push(
            format!("inventory.phases[{phase_index}].components[{component_index}].source_library"),
            ValidationIssueKind::InvalidValue,
            "source library must not be empty when provenance is provided",
        );
        None
    } else {
        Some(source_library.to_string())
    }
}

impl CandidatePolicyDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedCandidatePolicy> {
        let lower = parse_positive(
            &self.temperature_lower_k,
            "inventory.candidate_policy.temperature_lower_k",
            report,
        );
        let upper = parse_positive(
            &self.temperature_upper_k,
            "inventory.candidate_policy.temperature_upper_k",
            report,
        );
        if let (Some(lower), Some(upper)) = (lower, upper) {
            if lower > upper {
                report.push(
                    "inventory.candidate_policy.temperature_range",
                    ValidationIssueKind::InvalidValue,
                    "candidate temperature lower bound must not exceed upper bound",
                );
            }
        }
        let max_candidates = parse_usize_at_least(
            &self.max_candidates,
            "inventory.candidate_policy.max_candidates",
            1,
            report,
        );
        let states = self.physical_states.clone();
        if let (Some(lower), Some(upper), Some(max_candidates)) = (lower, upper, max_candidates) {
            if lower <= upper {
                Some(ValidatedCandidatePolicy {
                    element_mode: self.element_mode,
                    physical_states: states,
                    temperature_lower_k: lower,
                    temperature_upper_k: upper,
                    max_candidates,
                })
            } else {
                None
            }
        } else {
            None
        }
    }
}

impl EquilibriumLookupDraft {
    fn validate(&self, report: &mut EquilibriumGuiValidationReport) -> Option<ValidatedLookup> {
        match self {
            Self::Default => Some(ValidatedLookup::Default),
            Self::Explicit {
                priority_libraries,
                permitted_libraries,
                explicit_search_instructions,
                search_in_nist,
            } => {
                validate_library_names(priority_libraries, "lookup.priority_libraries", report);
                validate_library_names(permitted_libraries, "lookup.permitted_libraries", report);
                for (substance, library) in explicit_search_instructions {
                    if substance.trim().is_empty() || library.trim().is_empty() {
                        report.push(
                            "lookup.explicit_search_instructions",
                            ValidationIssueKind::InvalidValue,
                            "explicit lookup instruction keys and values must not be empty",
                        );
                    }
                }
                Some(ValidatedLookup::Explicit {
                    priority_libraries: priority_libraries.clone(),
                    permitted_libraries: permitted_libraries.clone(),
                    explicit_search_instructions: explicit_search_instructions.clone(),
                    search_in_nist: *search_in_nist,
                })
            }
        }
    }
}

fn validate_library_names(
    libraries: &[String],
    field: &str,
    report: &mut EquilibriumGuiValidationReport,
) {
    let mut seen = HashSet::new();
    for (index, library) in libraries.iter().enumerate() {
        let name = library.trim();
        if name.is_empty() {
            report.push(
                format!("{field}[{index}]"),
                ValidationIssueKind::MissingValue,
                "library name must not be empty",
            );
        } else if !seen.insert(name) {
            report.push(
                format!("{field}[{index}]"),
                ValidationIssueKind::Duplicate,
                format!("library '{name}' is listed more than once"),
            );
        }
    }
}

impl EquilibriumPhaseModeDraft {
    fn validate(&self, report: &mut EquilibriumGuiValidationReport) -> Option<ValidatedPhaseMode> {
        match self {
            Self::FixedDeclared => Some(ValidatedPhaseMode::FixedDeclared),
            Self::Bounded {
                phase_epsilon,
                dg_create,
                dg_keep,
                max_phase_iterations,
            } => {
                let epsilon = parse_positive(phase_epsilon, "phase_mode.phase_epsilon", report);
                let create = parse_finite(dg_create, "phase_mode.dg_create", report);
                let keep = parse_finite(dg_keep, "phase_mode.dg_keep", report);
                let iterations = parse_usize_at_least(
                    max_phase_iterations,
                    "phase_mode.max_phase_iterations",
                    1,
                    report,
                );
                match (epsilon, create, keep, iterations) {
                    (
                        Some(phase_epsilon),
                        Some(dg_create),
                        Some(dg_keep),
                        Some(max_phase_iterations),
                    ) => Some(ValidatedPhaseMode::Bounded {
                        phase_epsilon,
                        dg_create,
                        dg_keep,
                        max_phase_iterations,
                    }),
                    _ => None,
                }
            }
        }
    }
}

impl EquilibriumSolverConfigDraft {
    fn validate(&self, report: &mut EquilibriumGuiValidationReport) -> Option<ValidatedSolver> {
        if let EquilibriumSolverDraft::CustomCascade { backends } = &self.selection {
            if backends.is_empty() {
                report.push(
                    "solver.selection.backends",
                    ValidationIssueKind::MissingValue,
                    "custom solver cascade must contain at least one backend",
                );
            }
            let mut seen = HashSet::new();
            for (index, backend) in backends.iter().enumerate() {
                if !seen.insert(*backend) {
                    report.push(
                        format!("solver.selection.backends[{index}]"),
                        ValidationIssueKind::Duplicate,
                        format!(
                            "backend '{}' is repeated in the custom cascade",
                            backend.label()
                        ),
                    );
                }
            }
        }
        let tolerance = parse_optional_positive(
            &self.overrides.tolerance,
            "solver.overrides.tolerance",
            report,
        );
        let max_iterations = parse_optional_usize_at_least(
            &self.overrides.max_iterations,
            "solver.overrides.max_iterations",
            1,
            report,
        );
        let trace_seed_policy = self
            .overrides
            .trace_seed_policy
            .as_ref()
            .and_then(|policy| policy.validate(report));
        let cascade_budget = self
            .overrides
            .cascade_budget
            .as_ref()
            .and_then(|budget| budget.validate(report));
        if report.has_errors() && tolerance.is_none() && max_iterations.is_none() {
            // The caller still receives the complete report; this branch only
            // avoids manufacturing a validated override from malformed text.
        }
        Some(ValidatedSolver {
            selection: self.selection.clone(),
            tolerance,
            max_iterations,
            scaling_enabled: self.overrides.scaling_enabled,
            trace_seed_policy,
            cascade_budget,
        })
    }
}

impl GuiSolverCascadeBudgetDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedSolverCascadeBudget> {
        let max_attempts = parse_usize_at_least(
            &self.max_attempts,
            "solver.overrides.cascade_budget.max_attempts",
            1,
            report,
        );
        let max_iterations_per_attempt = parse_usize_at_least(
            &self.max_iterations_per_attempt,
            "solver.overrides.cascade_budget.max_iterations_per_attempt",
            1,
            report,
        );
        let max_total_iterations = parse_usize_at_least(
            &self.max_total_iterations,
            "solver.overrides.cascade_budget.max_total_iterations",
            1,
            report,
        );
        match (
            max_attempts,
            max_iterations_per_attempt,
            max_total_iterations,
        ) {
            (Some(max_attempts), Some(max_iterations_per_attempt), Some(max_total_iterations)) => {
                Some(ValidatedSolverCascadeBudget {
                    max_attempts,
                    max_iterations_per_attempt,
                    max_total_iterations,
                })
            }
            _ => None,
        }
    }
}

impl GuiTraceSeedPolicyDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedTraceSeedPolicy> {
        match self {
            Self::Absolute { floor } => {
                parse_positive(floor, "solver.overrides.trace_seed_policy.floor", report)
                    .map(|floor| ValidatedTraceSeedPolicy::Absolute { floor })
            }
            Self::RelativeToLargestInitialMole {
                fraction,
                minimum_floor,
            } => {
                let fraction = parse_positive(
                    fraction,
                    "solver.overrides.trace_seed_policy.fraction",
                    report,
                );
                let minimum_floor = parse_positive(
                    minimum_floor,
                    "solver.overrides.trace_seed_policy.minimum_floor",
                    report,
                );
                if let Some(value) = fraction {
                    if value > 1.0 {
                        report.push(
                            "solver.overrides.trace_seed_policy.fraction",
                            ValidationIssueKind::InvalidValue,
                            "fraction must not exceed 1",
                        );
                    }
                }
                match (fraction, minimum_floor) {
                    (Some(fraction), Some(minimum_floor)) if fraction <= 1.0 => {
                        Some(ValidatedTraceSeedPolicy::RelativeToLargestInitialMole {
                            fraction,
                            minimum_floor,
                        })
                    }
                    _ => None,
                }
            }
        }
    }
}

impl EquilibriumDiagnosticsDraft {
    fn validate(
        &self,
        _report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedDiagnostics> {
        Some(ValidatedDiagnostics {
            collect_timing: self.collect_timing,
            retain_backend_attempts: self.retain_backend_attempts,
            retain_conservation_report: self.retain_conservation_report,
            retain_phase_transitions: self.retain_phase_transitions,
            keq_validation: self.keq_validation,
        })
    }
}

impl EquilibriumPostprocessingDraft {
    fn validate(
        &self,
        report: &mut EquilibriumGuiValidationReport,
    ) -> Option<ValidatedPostprocessing> {
        let resampling = match &self.resampling {
            GuiResamplingDraft::None => None,
            GuiResamplingDraft::Pchip { output_points, .. } => parse_usize_at_least(
                output_points,
                "postprocessing.resampling.output_points",
                2,
                report,
            ),
        };
        Some(ValidatedPostprocessing {
            plot_target: self.plot_target,
            result_basis: self.result_basis,
            resampling,
            y_scale: self.y_scale,
        })
    }
}

fn phase_model_matches_state(state: GuiPhysicalState, model: GuiPhaseModel) -> bool {
    matches!(
        (state, model),
        (GuiPhysicalState::Gas, GuiPhaseModel::IdealGas)
            | (GuiPhysicalState::Liquid, GuiPhaseModel::PureCondensed)
            | (GuiPhysicalState::Solid, GuiPhaseModel::PureCondensed)
            | (GuiPhysicalState::Condensed, GuiPhaseModel::PureCondensed)
    )
}

fn parse_finite(
    raw: &str,
    field: &str,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<f64> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        report.push(
            field,
            ValidationIssueKind::MissingValue,
            "value is required",
        );
        return None;
    }
    match trimmed.parse::<f64>() {
        Ok(value) if value.is_finite() => Some(value),
        _ => {
            report.push(
                field,
                ValidationIssueKind::InvalidValue,
                "value must be a finite number",
            );
            None
        }
    }
}

fn parse_positive(
    raw: &str,
    field: &str,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<f64> {
    parse_finite(raw, field, report).and_then(|value| {
        if value > 0.0 {
            Some(value)
        } else {
            report.push(
                field,
                ValidationIssueKind::InvalidValue,
                "value must be greater than zero",
            );
            None
        }
    })
}

fn parse_non_negative(
    raw: &str,
    field: &str,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<f64> {
    parse_finite(raw, field, report).and_then(|value| {
        if value >= 0.0 {
            Some(value)
        } else {
            report.push(
                field,
                ValidationIssueKind::InvalidValue,
                "moles must be non-negative",
            );
            None
        }
    })
}

fn parse_usize_at_least(
    raw: &str,
    field: &str,
    minimum: usize,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<usize> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        report.push(
            field,
            ValidationIssueKind::MissingValue,
            "integer value is required",
        );
        return None;
    }
    match trimmed.parse::<usize>() {
        Ok(value) if value >= minimum => Some(value),
        Ok(_) => {
            report.push(
                field,
                ValidationIssueKind::InvalidValue,
                format!("value must be at least {minimum}"),
            );
            None
        }
        Err(_) => {
            report.push(
                field,
                ValidationIssueKind::InvalidValue,
                "value must be a positive integer",
            );
            None
        }
    }
}

fn parse_optional_positive(
    raw: &str,
    field: &str,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<f64> {
    if raw.trim().is_empty() {
        None
    } else {
        parse_positive(raw, field, report)
    }
}

fn parse_optional_usize_at_least(
    raw: &str,
    field: &str,
    minimum: usize,
    report: &mut EquilibriumGuiValidationReport,
) -> Option<usize> {
    if raw.trim().is_empty() {
        None
    } else {
        parse_usize_at_least(raw, field, minimum, report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagnostics_default_to_applicable_keq_validation_without_global_logging() {
        let diagnostics = EquilibriumDiagnosticsDraft::default();
        assert_eq!(
            diagnostics.keq_validation,
            GuiKeqValidationMode::WhenApplicable
        );
        assert!(!diagnostics.collect_timing);
    }

    #[test]
    fn cascade_budget_override_requires_positive_limits() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.solver.overrides.cascade_budget = Some(GuiSolverCascadeBudgetDraft {
            max_attempts: "0".into(),
            max_iterations_per_attempt: "10".into(),
            max_total_iterations: "20".into(),
        });
        let report = document
            .validate_for_run()
            .expect_err("zero cascade attempts must be rejected");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| { issue.field == "solver.overrides.cascade_budget.max_attempts" })
        );
    }

    #[test]
    fn default_document_is_runnable_and_roundtrips() {
        let document = EquilibriumGuiDocument::new();
        let validated = document
            .validate_for_run()
            .expect("default must be runnable");
        assert!(matches!(
            validated.problem,
            ValidatedProblem::FixedPt { .. }
        ));
        let json = document.to_json().expect("document serializes");
        let restored = EquilibriumGuiDocument::from_json(&json).expect("document restores");
        assert_eq!(document, restored);
    }

    #[test]
    fn simple_ideal_gas_preset_uses_the_canonical_multicomponent_model() {
        let document = EquilibriumGuiDocument::simple_ideal_gas(["H2", "O2", "N2"]);
        let EquilibriumInventoryDraft::ExplicitPhases { phases } = &document.config.inventory
        else {
            panic!("simple preset must use explicit phase inventory");
        };
        assert_eq!(phases.len(), 1);
        assert_eq!(phases[0].model, GuiPhaseModel::IdealGas);
        assert_eq!(phases[0].components.len(), 3);
        assert!(document.validate_for_run().is_ok());
    }

    #[test]
    fn range_validation_accepts_both_directions() {
        let mut document = EquilibriumGuiDocument::new();
        if let EquilibriumProblemDraft::FixedPt { temperature, .. } = &mut document.config.problem {
            *temperature = TemperatureDraft::Range {
                start_k: "300".into(),
                end_k: "1200".into(),
                point_count: "10".into(),
            };
        } else {
            panic!("default problem must be P,T");
        }
        assert!(document.validate_for_run().is_ok());
        if let EquilibriumProblemDraft::FixedPt { temperature, .. } = &mut document.config.problem {
            if let TemperatureDraft::Range { start_k, end_k, .. } = temperature {
                std::mem::swap(start_k, end_k);
            }
        } else {
            panic!("default problem must be P,T");
        }
        assert!(document.validate_for_run().is_ok());
    }

    #[test]
    fn invalid_numeric_fields_are_reported_by_field() {
        let mut document = EquilibriumGuiDocument::new();
        let EquilibriumProblemDraft::FixedPt {
            pressure_pa,
            temperature,
            ..
        } = &mut document.config.problem
        else {
            panic!("default problem must be P,T");
        };
        *pressure_pa = "NaN".into();
        *temperature = TemperatureDraft::Range {
            start_k: "1000".into(),
            end_k: "1000".into(),
            point_count: "1".into(),
        };
        let report = document
            .validate_for_run()
            .expect_err("invalid values must fail");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| issue.field == "problem.pressure_pa")
        );
        assert!(
            report
                .issues
                .iter()
                .any(|issue| issue.field == "temperature.point_count")
        );
        assert!(
            report
                .issues
                .iter()
                .any(|issue| issue.field == "temperature")
        );
    }

    #[test]
    fn trace_seed_override_validates_absolute_and_relative_policies() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.solver.overrides.trace_seed_policy =
            Some(GuiTraceSeedPolicyDraft::RelativeToLargestInitialMole {
                fraction: "1e-12".into(),
                minimum_floor: "1e-30".into(),
            });
        let validated = document
            .validate_for_run()
            .expect("valid relative trace policy must validate");
        assert!(matches!(
            validated.solver.trace_seed_policy,
            Some(ValidatedTraceSeedPolicy::RelativeToLargestInitialMole { .. })
        ));

        document.config.solver.overrides.trace_seed_policy =
            Some(GuiTraceSeedPolicyDraft::Absolute { floor: "0".into() });
        let report = document
            .validate_for_run()
            .expect_err("zero trace floor must fail");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| { issue.field == "solver.overrides.trace_seed_policy.floor" })
        );

        document.config.solver.overrides.trace_seed_policy =
            Some(GuiTraceSeedPolicyDraft::RelativeToLargestInitialMole {
                fraction: "2".into(),
                minimum_floor: "1e-30".into(),
            });
        let report = document
            .validate_for_run()
            .expect_err("relative fraction above one must fail");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| { issue.field == "solver.overrides.trace_seed_policy.fraction" })
        );
    }

    #[test]
    fn backend_catalog_is_exhaustive_and_has_unique_labels() {
        let labels = GuiSolverBackend::ALL
            .into_iter()
            .map(GuiSolverBackend::label)
            .collect::<Vec<_>>();
        assert_eq!(labels.len(), 9);
        let unique = labels.iter().collect::<HashSet<_>>();
        assert_eq!(unique.len(), labels.len());
        assert!(labels.iter().any(|label| label.contains("Legacy")));
        assert!(labels.iter().any(|label| label.contains("RST")));
    }

    #[test]
    fn same_substance_in_two_phases_is_valid_but_duplicate_in_one_phase_is_not() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.inventory = EquilibriumInventoryDraft::ExplicitPhases {
            phases: vec![
                PhaseDraft {
                    id: "gas".into(),
                    physical_state: GuiPhysicalState::Gas,
                    model: GuiPhaseModel::IdealGas,
                    components: vec![ComponentDraft {
                        substance: "H2O".into(),
                        initial_moles: "1".into(),
                        source_library: None,
                    }],
                },
                PhaseDraft {
                    id: "liquid".into(),
                    physical_state: GuiPhysicalState::Liquid,
                    model: GuiPhaseModel::PureCondensed,
                    components: vec![ComponentDraft {
                        substance: "H2O".into(),
                        initial_moles: "0.1".into(),
                        source_library: None,
                    }],
                },
            ],
        };
        assert!(document.validate_for_run().is_ok());

        if let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
        {
            let duplicate = phases[1].components[0].clone();
            phases[1].components.push(duplicate);
        }
        let report = document
            .validate_for_run()
            .expect_err("duplicate must fail");
        assert!(report.contains_kind(ValidationIssueKind::Duplicate));
    }

    /*
     * The remaining tests below intentionally stay after the structural tests
     * above. Keeping the document tests together makes schema changes easy to
     * review and prevents GUI code from becoming the only validation oracle.
     */
    #[test]
    fn fixed_ph_validates_and_roundtrips_as_total_enthalpy() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.problem = EquilibriumProblemDraft::FixedPh {
            pressure_pa: "101325".into(),
            reference_pressure_pa: "101325".into(),
            target_enthalpy_j: "10000".into(),
            temperature_bounds: PhTemperatureBoundsDraft::default(),
        };
        document
            .validate_for_run()
            .expect("valid P,H controls must reach request building");
        let json = document.to_json().expect("P,H document serializes");
        let restored = EquilibriumGuiDocument::from_json(&json).expect("P,H restores");
        assert_eq!(document, restored);
    }

    #[test]
    fn fixed_ph_rejects_a_seed_outside_its_temperature_bracket() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.problem = EquilibriumProblemDraft::FixedPh {
            pressure_pa: "101325".into(),
            reference_pressure_pa: "101325".into(),
            target_enthalpy_j: "10000".into(),
            temperature_bounds: PhTemperatureBoundsDraft {
                lower_k: "300".into(),
                upper_k: "1000".into(),
                seed_k: "1500".into(),
            },
        };
        let report = document
            .validate_for_run()
            .expect_err("seed outside the bracket must be rejected");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| issue.field == "problem.temperature_bounds.seed_k")
        );
    }

    #[test]
    fn element_mode_allows_preview_editing_but_requires_assignments_to_run() {
        let mut document = EquilibriumGuiDocument::new();
        document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
            elements: vec!["C".into(), "H".into(), "O".into()],
            candidate_policy: CandidatePolicyDraft::default(),
            assignments: Vec::new(),
        };
        assert!(document.validate().is_ok());
        let report = document
            .validate_for_run()
            .expect_err("unassigned candidates cannot run");
        assert!(
            report
                .issues
                .iter()
                .any(|issue| issue.field == "inventory.assignments")
        );
    }

    #[test]
    fn incompatible_phase_model_is_rejected() {
        let mut document = EquilibriumGuiDocument::new();
        let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut document.config.inventory
        else {
            panic!("default inventory must be explicit");
        };
        phases[0].physical_state = GuiPhysicalState::Solid;
        let report = document
            .validate_for_run()
            .expect_err("solid ideal gas is invalid");
        assert!(report.contains_kind(ValidationIssueKind::Incompatible));
    }

    #[test]
    fn future_schema_is_rejected_before_deserialization_can_run() {
        let mut document = EquilibriumGuiDocument::new();
        document.schema_version += 1;
        let json = serde_json::to_string(&document).expect("document serializes");
        let error = EquilibriumGuiDocument::from_json(&json).expect_err("future schema must fail");
        assert!(matches!(
            error,
            EquilibriumGuiDocumentError::InvalidSchema(_)
        ));
    }
}
