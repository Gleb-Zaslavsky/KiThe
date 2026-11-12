//! Reactor IVP document parsing.
//!
//! This module owns the reactor physics layer for condensed-phase IVP tasks.
//! Solver-specific settings are delegated to RustedSciThe's IVP parser, so this
//! file only parses the physical reactor contract, the typed initial state, and
//! then bundles the solver settings as a typed RustedSciThe spec.

use crate::Kinetics::mechfinder_api::{ReactionData, ReactionKinetics};
use crate::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use crate::SimpleReactorsIVP::SimpleReactorIVP::{IvpError, ReactorIvpPhysicalConfig};
use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
use RustedSciThe::command_interpreter::task_parser_ivp::{
    IvpSolverSettingsSpec, Lsode2TaskExecutionSpec, parse_ivp_solver_settings_from_document,
};
use RustedSciThe::command_interpreter::task_parser_ivp_hybrid::{
    HybridStiffSettings, parse_hybrid_stiff_settings_from_document,
    serialize_hybrid_solver_fragment,
};
use RustedSciThe::numerical::LSODE2::config::{
    Lsode2AotProfile, Lsode2AotToolchain, Lsode2LinearSolverChoice, Lsode2LinearSolverPolicy,
    Lsode2LinearSystemStructure, Lsode2NativeExecutionConfig,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

type GenericSectionMap = HashMap<String, Option<Vec<Value>>>;

/// Normalized reactor IVP document contract.
///
/// The physics and local reactor-only controls stay in KiThe, while solver
/// settings remain in the typed RustedSciThe spec that parsed them.
#[derive(Debug, Clone)]
pub struct ReactorIvpDocumentSpec {
    pub problem_name: Option<String>,
    pub problem_description: Option<String>,
    pub physical: ReactorIvpPhysicalConfig,
    pub scaling: ScalingConfig,
    /// Lower dimensionless spatial coordinate z=x/L (legacy document key `x0`).
    pub x0: f64,
    /// Upper dimensionless spatial coordinate z=x/L (legacy document key `x_bound`).
    pub x_bound: f64,
    pub initial_conditions: HashMap<String, f64>,
    pub thermal_effects: Vec<f64>,
    pub mechanism: Option<ReactorIvpMechanismSpec>,
    pub stop_condition: ReactorIvpStopConditionSpec,
    pub solver_settings: IvpSolverSettingsSpec,
    pub hybrid_solver_settings: Option<HybridStiffSettings>,
    pub hybrid_solver_fragment: Option<String>,
}

/// Optional editable reaction mechanism payload stored with an IVP task.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactorIvpMechanismSpec {
    pub substances: Vec<String>,
    pub reactions: Vec<ReactionData>,
}

/// Local stop-condition contract for the condensed-phase IVP task document.
///
/// The solver backend already knows how to execute the stop rule, but the
/// document parser keeps the user-facing selection here because it belongs to
/// the reactor task rather than the generic RustedSciThe solver syntax.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorIvpStopConditionSpec {
    pub enabled: bool,
    pub species: Option<String>,
    pub threshold: Option<f64>,
}

impl Default for ReactorIvpStopConditionSpec {
    fn default() -> Self {
        Self {
            enabled: false,
            species: None,
            threshold: None,
        }
    }
}

impl ReactorIvpStopConditionSpec {
    /// Validate the stop-condition contract against a known species order.
    pub fn validate_against_species(&self, species_order: &[String]) -> Result<(), IvpError> {
        if !self.enabled {
            return Ok(());
        }

        let species = self.species.as_ref().ok_or_else(|| {
            IvpError::MissingData("Missing `species` in `stop_condition` section".to_string())
        })?;
        if !species_order.iter().any(|name| name == species) {
            return Err(IvpError::InvalidConfiguration(format!(
                "Unknown stop-condition species `{species}`"
            )));
        }

        let threshold = self.threshold.ok_or_else(|| {
            IvpError::MissingData("Missing `threshold` in `stop_condition` section".to_string())
        })?;
        if !threshold.is_finite() || threshold < 0.0 {
            return Err(IvpError::InvalidNumericValue(format!(
                "stop_condition.threshold must be finite and non-negative, got {}",
                threshold
            )));
        }

        Ok(())
    }
}

impl ReactorIvpDocumentSpec {
    /// Validate the full parsed document without mutating any task state.
    ///
    /// This is the pure validation step used before task application or
    /// equation preview. It keeps parse/validate/apply as separate phases.
    pub fn validate(&self, species_order: &[String]) -> Result<(), IvpError> {
        self.validate_against_species(species_order)
    }

    /// Validate the initial-state entries against a known substance order.
    pub fn validate_against_species(&self, species_order: &[String]) -> Result<(), IvpError> {
        let effective_species = self
            .mechanism
            .as_ref()
            .map(|mechanism| mechanism.substances.as_slice())
            .unwrap_or(species_order);
        if let Some(mechanism) = &self.mechanism {
            let mut unique = std::collections::HashSet::new();
            for substance in &mechanism.substances {
                if substance.trim().is_empty() || !unique.insert(substance) {
                    return Err(IvpError::InvalidConfiguration(format!(
                        "mechanism contains an empty or duplicate substance `{substance}`"
                    )));
                }
            }
            if mechanism.reactions.len() != self.thermal_effects.len() {
                return Err(IvpError::InvalidConfiguration(format!(
                    "mechanism reaction count ({}) must match thermal_effects count ({})",
                    mechanism.reactions.len(),
                    self.thermal_effects.len()
                )));
            }
            for reaction in &mechanism.reactions {
                reaction.validate().map_err(|error| {
                    IvpError::InvalidConfiguration(format!(
                        "invalid reaction `{}`: {error}",
                        reaction.eq
                    ))
                })?;
            }
        }
        for required in effective_species {
            if !self.initial_conditions.contains_key(required) {
                return Err(IvpError::MissingData(format!(
                    "Missing initial condition for substance `{required}`"
                )));
            }
        }
        if !self.initial_conditions.contains_key("T") {
            return Err(IvpError::MissingData(
                "Missing `T` in initial conditions".to_string(),
            ));
        }
        if !self.initial_conditions.contains_key("q") {
            return Err(IvpError::MissingData(
                "Missing `q` in initial conditions".to_string(),
            ));
        }
        self.stop_condition
            .validate_against_species(effective_species)?;
        Ok(())
    }

    /// Apply the parsed physics contract to a reactor task.
    ///
    /// The solver settings remain in the returned document spec so the caller
    /// can decide when to translate them into a concrete backend config.
    pub fn apply_to_task(
        &self,
        task: &mut crate::SimpleReactorsIVP::SimpleReactorIVP::SimpleReactorTask,
    ) -> Result<(), IvpError> {
        if let Some(mechanism) = &self.mechanism {
            let groups = task.kindata.groups().cloned();
            task.kindata
                .set_reaction_data_directly(mechanism.reactions.clone(), groups)
                .map_err(|error| {
                    IvpError::InvalidConfiguration(format!(
                        "failed to install IVP reaction mechanism: {error}"
                    ))
                })?;
            task.kindata.substances = mechanism.substances.clone();
            task.kindata
                .stecheodata
                .set_substances(mechanism.substances.clone());
            if !mechanism.reactions.is_empty() {
                task.kindata.analyze_reactions().map_err(|error| {
                    IvpError::InvalidConfiguration(format!(
                        "failed to analyze IVP reaction mechanism: {error}"
                    ))
                })?;
            }
        }
        if let Some(name) = &self.problem_name {
            task.set_problem_name(name);
        }
        if let Some(description) = &self.problem_description {
            task.set_problem_description(description);
        }
        task.set_density(self.physical.ro)?;
        task.set_transport_properties(self.physical.Lambda, self.physical.Cp);
        task.m = self.physical.m;
        task.L = self.physical.L;
        task.set_scaling(self.scaling.clone())?;
        task.set_initial_conditions(self.initial_conditions.clone());
        task.set_thermal_effects(self.thermal_effects.clone());
        task.set_hybrid_solver_settings(self.hybrid_solver_settings.clone());
        Ok(())
    }

    /// Serialize the normalized reactor-IVP document back into parser-ready text.
    ///
    /// The output is canonical enough for load-save-load round trips even when
    /// optional solver settings were originally supplied in a shorter form.
    pub fn document_to_string(&self) -> String {
        let mut out = String::new();

        let mut process_conditions = Vec::new();
        if let Some(name) = &self.problem_name {
            process_conditions.push(("problem_name", name.clone()));
        }
        if let Some(description) = &self.problem_description {
            process_conditions.push(("problem_description", description.clone()));
        }
        process_conditions.push(("ro", self.physical.ro.to_string()));
        process_conditions.push(("Cp", self.physical.Cp.to_string()));
        process_conditions.push(("Lambda", self.physical.Lambda.to_string()));
        process_conditions.push(("m", self.physical.m.to_string()));
        process_conditions.push(("L", self.physical.L.to_string()));
        process_conditions.push(("thermal_effects", format_vector_f64(&self.thermal_effects)));
        push_section(&mut out, "process_conditions", &process_conditions);

        if let Some(mechanism) = &self.mechanism {
            let mut kinetics = vec![("substances".to_string(), mechanism.substances.join(", "))];
            for (index, reaction) in mechanism.reactions.iter().enumerate() {
                let ReactionKinetics::Elementary(data) = &reaction.data else {
                    continue;
                };
                kinetics.push((format!("reaction_{index:03}_equation"), reaction.eq.clone()));
                kinetics.push((
                    format!("reaction_{index:03}_arrhenius"),
                    format_vector_f64(&data.Arrenius),
                ));
            }
            let kinetics = kinetics
                .iter()
                .map(|(key, value)| (key.as_str(), value.clone()))
                .collect::<Vec<_>>();
            push_section(&mut out, "kinetics", &kinetics);
        }

        push_section(
            &mut out,
            "scaling",
            &[
                ("dT", self.scaling.dT.to_string()),
                ("T_scale", self.scaling.T_scale.to_string()),
                ("L", self.scaling.L.to_string()),
            ],
        );

        out.push_str("# Legacy x0/x_bound keys contain dimensionless z=x/L coordinates\n");
        push_section(
            &mut out,
            "integration_domain",
            &[
                ("x0", self.x0.to_string()),
                ("x_bound", self.x_bound.to_string()),
            ],
        );

        let mut initial_conditions = vec![
            (
                "T".to_string(),
                self.initial_conditions
                    .get("T")
                    .copied()
                    .unwrap_or(0.0)
                    .to_string(),
            ),
            (
                "q".to_string(),
                self.initial_conditions
                    .get("q")
                    .copied()
                    .unwrap_or(0.0)
                    .to_string(),
            ),
        ];
        let mut species_entries = self
            .initial_conditions
            .iter()
            .filter(|(name, _)| name.as_str() != "T" && name.as_str() != "q")
            .map(|(name, value)| (name.clone(), *value))
            .collect::<Vec<_>>();
        species_entries.sort_by(|left, right| left.0.cmp(&right.0));
        for (name, value) in species_entries {
            initial_conditions.push((name, value.to_string()));
        }
        push_section(&mut out, "initial_conditions", &initial_conditions);

        let mut stop_condition = vec![("enabled", self.stop_condition.enabled.to_string())];
        if let Some(species) = &self.stop_condition.species {
            stop_condition.push(("species", species.clone()));
        }
        if let Some(threshold) = self.stop_condition.threshold {
            stop_condition.push(("threshold", threshold.to_string()));
        }
        push_section(&mut out, "stop_condition", &stop_condition);

        if let Some(fragment) = &self.hybrid_solver_fragment {
            out.push_str(fragment);
            return out;
        }

        push_section(
            &mut out,
            "task",
            &[
                ("solver", "IVP".to_string()),
                (
                    "method",
                    method_document_name(&self.solver_settings.solver.method),
                ),
            ],
        );

        let mut solver_options = Vec::new();
        if let Some(value) = self.solver_settings.solver_options.max_iterations {
            solver_options.push(("max_iterations", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.rtol {
            solver_options.push(("rtol", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.atol {
            solver_options.push(("atol", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.max_step {
            solver_options.push(("max_step", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.first_step {
            solver_options.push(("first_step", format!("Some({value})")));
        }
        if let Some(value) = self.solver_settings.solver_options.vectorized {
            solver_options.push(("vectorized", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.parallel {
            solver_options.push(("parallel", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.neighborhood_check {
            solver_options.push(("neighborhood_check", value.to_string()));
        }
        if let Some(value) = self.solver_settings.solver_options.telemetry {
            let label = match value {
                RustedSciThe::command_interpreter::task_parser_ivp::IvpTelemetryModeSpec::Off => "off",
                RustedSciThe::command_interpreter::task_parser_ivp::IvpTelemetryModeSpec::Counters => "counters",
                RustedSciThe::command_interpreter::task_parser_ivp::IvpTelemetryModeSpec::Timings => "timings",
                RustedSciThe::command_interpreter::task_parser_ivp::IvpTelemetryModeSpec::Detailed => "detailed",
            };
            solver_options.push(("telemetry", label.to_string()));
        }
        if let Some(level) = self.solver_settings.solver_options.logging {
            solver_options.push(("logging", format!("{:?}", level).to_ascii_lowercase()));
        }

        if let Some(lsode2) = self.solver_settings.solver_options.lsode2.as_ref() {
            if let Some(controller) = lsode2.controller {
                solver_options.push(("lsode2_method_family", controller.mode.label().to_string()));
            }
            if let Some(assembly) = lsode2.symbolic_assembly {
                solver_options.push(("lsode2_symbolic_assembly", assembly.label().to_string()));
            }
            if let Some(execution) = &lsode2.symbolic_execution {
                match execution {
                    Lsode2TaskExecutionSpec::LambdifyExpr => {
                        solver_options
                            .push(("lsode2_symbolic_execution", "LambdifyExpr".to_string()));
                    }
                    Lsode2TaskExecutionSpec::Aot {
                        toolchain,
                        profile,
                        output_parent_dir,
                    } => {
                        solver_options.push(("lsode2_symbolic_execution", "AOT".to_string()));
                        solver_options.push((
                            "lsode2_aot_toolchain",
                            aot_toolchain_label(*toolchain).to_string(),
                        ));
                        solver_options.push((
                            "lsode2_aot_profile",
                            aot_profile_label(*profile).to_string(),
                        ));
                        if let Some(path) = output_parent_dir {
                            solver_options
                                .push(("lsode2_aot_output_dir", path.display().to_string()));
                        }
                    }
                }
            }
            if let Some(structure) = lsode2.linear_system_structure {
                match structure {
                    Lsode2LinearSystemStructure::Dense => {
                        solver_options.push(("lsode2_linear_structure", "dense".to_string()));
                    }
                    Lsode2LinearSystemStructure::Sparse => {
                        solver_options.push(("lsode2_linear_structure", "sparse".to_string()));
                    }
                    Lsode2LinearSystemStructure::Banded { kl, ku } => {
                        solver_options.push(("lsode2_linear_structure", "banded".to_string()));
                        solver_options.push(("lsode2_banded_kl", kl.to_string()));
                        solver_options.push(("lsode2_banded_ku", ku.to_string()));
                    }
                    Lsode2LinearSystemStructure::BandedAuto => {
                        solver_options.push(("lsode2_linear_structure", "banded".to_string()));
                    }
                }
            }
            if let Some(policy) = lsode2.linear_solver_policy {
                match policy {
                    Lsode2LinearSolverPolicy::Auto => {
                        solver_options.push(("lsode2_linear_solver_policy", "auto".to_string()));
                    }
                    Lsode2LinearSolverPolicy::Force(choice) => {
                        solver_options.push((
                            "lsode2_linear_solver_policy",
                            linear_solver_choice_label(choice).to_string(),
                        ));
                    }
                }
            }
            if let Some(native_execution) = lsode2.native_execution {
                match native_execution {
                    Lsode2NativeExecutionConfig::BridgeSolve => {
                        solver_options
                            .push(("lsode2_native_execution", "bridge_solve".to_string()));
                    }
                    Lsode2NativeExecutionConfig::Disabled => {
                        solver_options
                            .push(("lsode2_native_execution", "bridge_solve".to_string()));
                    }
                    Lsode2NativeExecutionConfig::ProbeBeforeBridge {
                        max_step_attempts,
                        max_accepted_steps,
                    } => {
                        solver_options
                            .push(("lsode2_native_execution", "probe_before_bridge".to_string()));
                        solver_options.push((
                            "lsode2_native_max_step_attempts",
                            max_step_attempts.to_string(),
                        ));
                        solver_options.push((
                            "lsode2_native_max_accepted_steps",
                            max_accepted_steps.to_string(),
                        ));
                    }
                    Lsode2NativeExecutionConfig::NativeSolve {
                        max_step_attempts,
                        max_accepted_steps,
                    } => {
                        solver_options
                            .push(("lsode2_native_execution", "faithful_bdf_solve".to_string()));
                        solver_options.push((
                            "lsode2_native_max_step_attempts",
                            max_step_attempts.to_string(),
                        ));
                        solver_options.push((
                            "lsode2_native_max_accepted_steps",
                            max_accepted_steps.to_string(),
                        ));
                    }
                }
            }
        }

        if solver_options.is_empty()
            && self.solver_settings.solver.method
                == RustedSciThe::command_interpreter::task_parser_ivp::IvpMethodSpec::Lsode2
        {
            solver_options.push(("lsode2_method_family", "auto".to_string()));
            solver_options.push(("lsode2_symbolic_assembly", "AtomView".to_string()));
            solver_options.push(("lsode2_symbolic_execution", "LambdifyExpr".to_string()));
            solver_options.push(("lsode2_linear_structure", "sparse".to_string()));
            solver_options.push(("lsode2_linear_solver_policy", "auto".to_string()));
            solver_options.push(("lsode2_native_execution", "bridge_solve".to_string()));
        }
        push_section(&mut out, "solver_options", &solver_options);

        out
    }
}

/// Parse a reactor IVP document from a normalized document map.
pub fn parse_reactor_ivp_document(
    document: &DocumentMap,
) -> Result<ReactorIvpDocumentSpec, IvpError> {
    let process_conditions = required_section(document, "process_conditions")?;
    let physical = ReactorIvpPhysicalConfig {
        ro: required_f64(process_conditions, "process_conditions", "ro")?,
        Cp: required_f64(process_conditions, "process_conditions", "Cp")?,
        Lambda: required_f64(process_conditions, "process_conditions", "Lambda")?,
        m: required_f64(process_conditions, "process_conditions", "m")?,
        L: required_f64(process_conditions, "process_conditions", "L")?,
    };

    let problem_name = optional_string(process_conditions, "process_conditions", "problem_name")?;
    let problem_description = optional_string(
        process_conditions,
        "process_conditions",
        "problem_description",
    )?;
    let scaling = parse_scaling(document, physical.L)?;
    let (x0, x_bound) = parse_integration_domain(document)?;
    let initial_conditions = parse_initial_conditions(document)?;
    let thermal_effects =
        required_vector_f64(process_conditions, "process_conditions", "thermal_effects")?;
    let mechanism = parse_mechanism(document)?;
    let stop_condition = parse_stop_condition(document)?;
    let solver_settings = parse_ivp_solver_settings_from_document(document)
        .map_err(|error| IvpError::InvalidConfiguration(error.to_string()))?;
    let hybrid_solver_settings = parse_hybrid_stiff_settings_from_document(document)
        .map_err(|error| IvpError::InvalidConfiguration(error.to_string()))?;
    let hybrid_solver_fragment = serialize_hybrid_solver_fragment(document)
        .map_err(|error| IvpError::InvalidConfiguration(error.to_string()))?;
    if stop_condition.enabled
        && solver_settings.solver.method
            != RustedSciThe::command_interpreter::task_parser_ivp::IvpMethodSpec::Lsode2
    {
        return Err(IvpError::InvalidConfiguration(
            "the condensed species stop condition is currently supported only by LSODE2".into(),
        ));
    }

    Ok(ReactorIvpDocumentSpec {
        problem_name,
        problem_description,
        physical,
        scaling,
        x0,
        x_bound,
        initial_conditions,
        thermal_effects,
        mechanism,
        stop_condition,
        solver_settings,
        hybrid_solver_settings: Some(hybrid_solver_settings),
        hybrid_solver_fragment: Some(hybrid_solver_fragment),
    })
}

fn method_document_name(
    method: &RustedSciThe::command_interpreter::task_parser_ivp::IvpMethodSpec,
) -> String {
    use RustedSciThe::command_interpreter::task_parser_ivp::IvpMethodSpec;
    match method {
        IvpMethodSpec::Lsode2 => "LSODE2".to_string(),
        IvpMethodSpec::Bdf => "BDF".to_string(),
        IvpMethodSpec::Radau5 => "Radau5".to_string(),
        IvpMethodSpec::BackwardEuler => "BackwardEuler".to_string(),
        IvpMethodSpec::Lsode => "LSODE".to_string(),
        IvpMethodSpec::Lsoda => "LSODA".to_string(),
        IvpMethodSpec::NonStiff(name) => name.clone(),
    }
}

/// Parse a reactor IVP document directly from text.
pub fn parse_reactor_ivp_document_from_str(
    input: &str,
) -> Result<ReactorIvpDocumentSpec, IvpError> {
    let mut parser = DocumentParser::new(input.to_string());
    parser.parse_document().map_err(|error| {
        IvpError::InvalidConfiguration(format!("failed to parse reactor IVP document: {error}"))
    })?;
    let document = parser.get_result().ok_or_else(|| {
        IvpError::MissingData("reactor IVP document map was not produced".to_string())
    })?;
    parse_reactor_ivp_document(document)
}

/// Build a parser-ready canonical condensed-reactor IVP task template.
///
/// Physical inputs stay local to KiThe, while solver options use the typed
/// vocabulary currently supported by RustedSciThe for the selected family.
pub fn build_canonical_condensed_reactor_ivp_task_template() -> String {
    let mut template = String::new();

    // Task shell.
    template.push_str("task\n");
    template.push_str("solver: IVP\n");
    template.push_str("method: LSODE2\n\n");

    // Reactor physics stays in KiThe.
    template.push_str("process_conditions\n");
    template.push_str("problem_name: CondensedBurn\n");
    template.push_str("problem_description: Canonical condensed-phase IVP template\n");
    template.push_str("ro: 1800.0\n");
    template.push_str("Cp: 1200.0\n");
    template.push_str("Lambda: 0.35\n");
    template.push_str("# Mass flux per unit area [kg/(m^2*s)]\n");
    template.push_str("m: 0.02\n");
    template.push_str("L: 0.02\n");
    template.push_str("# Molar heat of reaction for each reaction [J/mol]\n");
    template.push_str("thermal_effects: [125000.0]\n\n");

    template.push_str("scaling\n");
    template.push_str("dT: 100.0\n");
    template.push_str("T_scale: 75.0\n");
    template.push_str("L: 0.02\n\n");

    template.push_str("integration_domain\n");
    template.push_str("# Legacy x0/x_bound keys contain dimensionless z=x/L coordinates\n");
    template.push_str("x0: 0.0\n");
    template.push_str("x_bound: 1.0\n\n");

    template.push_str("initial_conditions\n");
    template.push_str("T: 900.0\n");
    template.push_str("# Physical heat flux [W/m^2]; scaled internally by L/T_scale\n");
    template.push_str("q: 0.0\n");
    template.push_str("A: 1.0\n");
    template.push_str("B: 0.0\n\n");

    template.push_str("stop_condition\n");
    template.push_str("enabled: false\n\n");

    // Solver settings are deliberately parser-ready for RustedSciThe.
    template.push_str("solver_options\n");
    template.push_str("rtol: 1e-6\n");
    template.push_str("atol: 1e-8\n");
    template.push_str("max_step: 1e-3\n");
    template.push_str("first_step: Some(1e-6)\n");
    template.push_str("lsode2_method_family: auto\n");
    template.push_str("lsode2_symbolic_assembly: AtomView\n");
    template.push_str("lsode2_symbolic_execution: LambdifyExpr\n");
    template.push_str("lsode2_linear_structure: sparse\n");
    template.push_str("lsode2_linear_solver_policy: auto\n");

    template
}

fn parse_scaling(
    document: &DocumentMap,
    characteristic_length: f64,
) -> Result<ScalingConfig, IvpError> {
    let Some(section) = document.get("scaling") else {
        return Ok(ScalingConfig::new(100.0, characteristic_length, 100.0));
    };

    let dT = optional_f64(section, "scaling", "dT")?.unwrap_or(100.0);
    let T_scale = optional_f64(section, "scaling", "T_scale")?.unwrap_or(100.0);
    let L = optional_f64(section, "scaling", "L")?.unwrap_or(characteristic_length);

    if (L - characteristic_length).abs() > f64::EPSILON {
        return Err(IvpError::InvalidConfiguration(format!(
            "scaling.L ({L}) must match process_conditions.L ({characteristic_length})"
        )));
    }

    let scaling = ScalingConfig::new(dT, L, T_scale);
    scaling.validate().map_err(IvpError::from)?;
    Ok(scaling)
}

fn parse_integration_domain(document: &DocumentMap) -> Result<(f64, f64), IvpError> {
    let Some(section) = document.get("integration_domain") else {
        return Ok((0.0, 1.0));
    };

    let x0 = optional_f64(section, "integration_domain", "x0")?.unwrap_or(0.0);
    let x_bound = optional_f64(section, "integration_domain", "x_bound")?.unwrap_or(1.0);

    if !x0.is_finite() || !x_bound.is_finite() {
        return Err(IvpError::InvalidNumericValue(
            "integration_domain values must be finite".to_string(),
        ));
    }
    if x_bound <= x0 {
        return Err(IvpError::InvalidConfiguration(format!(
            "integration_domain.x_bound must be larger than x0 (x0={}, x_bound={})",
            x0, x_bound
        )));
    }

    Ok((x0, x_bound))
}

fn parse_initial_conditions(document: &DocumentMap) -> Result<HashMap<String, f64>, IvpError> {
    let section = required_section(document, "initial_conditions")?;
    let mut map = HashMap::with_capacity(section.len());
    for (field, values) in section {
        let values = values
            .as_ref()
            .ok_or_else(|| missing_field("initial_conditions", field))?;
        let value = values
            .first()
            .ok_or_else(|| missing_field("initial_conditions", field))?;
        let scalar =
            value_to_f64(value).ok_or_else(|| missing_field("initial_conditions", field))?;
        if !scalar.is_finite() {
            return Err(invalid_numeric("initial_conditions", field, scalar));
        }
        map.insert(field.clone(), scalar);
    }
    if !map.contains_key("T") {
        return Err(IvpError::MissingData(
            "Missing `T` in initial_conditions".to_string(),
        ));
    }
    if !map.contains_key("q") {
        return Err(IvpError::MissingData(
            "Missing `q` in initial_conditions".to_string(),
        ));
    }
    Ok(map)
}

fn parse_stop_condition(document: &DocumentMap) -> Result<ReactorIvpStopConditionSpec, IvpError> {
    let Some(section) = document.get("stop_condition") else {
        return Ok(ReactorIvpStopConditionSpec::default());
    };

    let enabled = optional_bool(section, "stop_condition", "enabled")?
        .unwrap_or_else(|| section.contains_key("species") || section.contains_key("threshold"));
    let species = optional_string(section, "stop_condition", "species")?;
    let threshold = optional_f64(section, "stop_condition", "threshold")?;

    Ok(ReactorIvpStopConditionSpec {
        enabled,
        species,
        threshold,
    })
}

fn parse_mechanism(document: &DocumentMap) -> Result<Option<ReactorIvpMechanismSpec>, IvpError> {
    let Some(section) = document.get("kinetics") else {
        return Ok(None);
    };

    if section.contains_key("mechanism_json") {
        let encoded = required_string(section, "kinetics", "mechanism_json")?;
        let json = encoded
            .trim()
            .strip_prefix('(')
            .and_then(|value| value.strip_suffix(')'))
            .ok_or_else(|| {
                IvpError::InvalidConfiguration(
                    "kinetics.mechanism_json must be wrapped in parentheses".to_string(),
                )
            })?;
        let mechanism = serde_json::from_str(json).map_err(|error| {
            IvpError::InvalidConfiguration(format!(
                "invalid kinetics.mechanism_json payload: {error}"
            ))
        })?;
        return Ok(Some(mechanism));
    }

    let Some(substance_values) = section.get("substances") else {
        return Ok(None);
    };
    let substance_values = substance_values
        .as_ref()
        .ok_or_else(|| missing_field("kinetics", "substances"))?;
    let substances = substance_values
        .iter()
        .map(|value| {
            value
                .as_string()
                .cloned()
                .ok_or_else(|| missing_field("kinetics", "substances"))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut reaction_indices = section
        .keys()
        .filter_map(|key| {
            key.strip_prefix("reaction_")?
                .strip_suffix("_equation")?
                .parse::<usize>()
                .ok()
        })
        .collect::<Vec<_>>();
    reaction_indices.sort_unstable();
    reaction_indices.dedup();

    let mut reactions = Vec::with_capacity(reaction_indices.len());
    for (expected_index, index) in reaction_indices.into_iter().enumerate() {
        if index != expected_index {
            return Err(IvpError::InvalidConfiguration(format!(
                "kinetics reaction indices must be contiguous from 000; found {index:03} at position {expected_index}"
            )));
        }
        let equation = required_string(
            section,
            "kinetics",
            &format!("reaction_{index:03}_equation"),
        )?;
        let arrhenius = required_vector_f64(
            section,
            "kinetics",
            &format!("reaction_{index:03}_arrhenius"),
        )?;
        if arrhenius.len() != 3 {
            return Err(IvpError::InvalidConfiguration(format!(
                "kinetics reaction {index:03} requires exactly A, n, and E Arrhenius values"
            )));
        }
        reactions.push(
            ReactionData::try_new_elementary(equation.clone(), arrhenius, None).map_err(
                |error| {
                    IvpError::InvalidConfiguration(format!(
                        "invalid kinetics reaction {index:03} `{equation}`: {error}"
                    ))
                },
            )?,
        );
    }

    if substances.is_empty() && reactions.is_empty() {
        return Ok(None);
    }
    if substances.is_empty() || reactions.is_empty() {
        return Err(IvpError::InvalidConfiguration(
            "readable kinetics requires both substances and at least one reaction".to_string(),
        ));
    }

    Ok(Some(ReactorIvpMechanismSpec {
        substances,
        reactions,
    }))
}

fn required_section<'a>(
    document: &'a DocumentMap,
    section: &str,
) -> Result<&'a GenericSectionMap, IvpError> {
    document
        .get(section)
        .ok_or_else(|| IvpError::MissingData(format!("Missing `{section}` section")))
}

fn required_f64(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<f64, IvpError> {
    let value = required_value(section, section_name, field)?;
    let scalar = value_to_f64(value).ok_or_else(|| missing_field(section_name, field))?;
    if !scalar.is_finite() {
        return Err(invalid_numeric(section_name, field, scalar));
    }
    Ok(scalar)
}

fn required_string(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<String, IvpError> {
    required_value(section, section_name, field)?
        .as_string()
        .cloned()
        .ok_or_else(|| missing_field(section_name, field))
}

fn optional_f64(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<Option<f64>, IvpError> {
    let Some(values) = section.get(field) else {
        return Ok(None);
    };
    let values = values
        .as_ref()
        .ok_or_else(|| missing_field(section_name, field))?;
    let Some(value) = values.first() else {
        return Ok(None);
    };
    let scalar = value_to_f64(value).ok_or_else(|| missing_field(section_name, field))?;
    if !scalar.is_finite() {
        return Err(invalid_numeric(section_name, field, scalar));
    }
    Ok(Some(scalar))
}

fn optional_string(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<Option<String>, IvpError> {
    let Some(values) = section.get(field) else {
        return Ok(None);
    };
    let values = values
        .as_ref()
        .ok_or_else(|| missing_field(section_name, field))?;
    let Some(value) = values.first() else {
        return Ok(None);
    };
    Ok(value
        .as_option_string()
        .cloned()
        .or_else(|| value.as_string().cloned()))
}

fn optional_bool(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<Option<bool>, IvpError> {
    let Some(values) = section.get(field) else {
        return Ok(None);
    };
    let values = values
        .as_ref()
        .ok_or_else(|| missing_field(section_name, field))?;
    let Some(value) = values.first() else {
        return Ok(None);
    };
    value
        .as_boolean()
        .ok_or_else(|| missing_field(section_name, field))
        .map(Some)
}

fn required_vector_f64(
    section: &GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<Vec<f64>, IvpError> {
    let value = required_value(section, section_name, field)?;
    let vector = value
        .as_vector()
        .ok_or_else(|| missing_field(section_name, field))?;
    if let Some(invalid) = vector.iter().copied().find(|value| !value.is_finite()) {
        return Err(invalid_numeric(section_name, field, invalid));
    }
    Ok(vector.clone())
}

fn required_value<'a>(
    section: &'a GenericSectionMap,
    section_name: &str,
    field: &str,
) -> Result<&'a Value, IvpError> {
    section
        .get(field)
        .ok_or_else(|| missing_field(section_name, field))?
        .as_ref()
        .ok_or_else(|| missing_field(section_name, field))?
        .first()
        .ok_or_else(|| missing_field(section_name, field))
}

fn value_to_f64(value: &Value) -> Option<f64> {
    value
        .as_float()
        .or_else(|| value.as_integer().map(|v| v as f64))
}

fn missing_field(section: &str, field: &str) -> IvpError {
    IvpError::MissingData(format!("Missing `{field}` in `{section}` section"))
}

fn invalid_numeric(section: &str, field: &str, value: f64) -> IvpError {
    IvpError::InvalidNumericValue(format!("`{section}.{field} = {value}` is not finite"))
}

fn push_section<K: AsRef<str>>(out: &mut String, section: &str, fields: &[(K, String)]) {
    if fields.is_empty() {
        return;
    }
    out.push_str(section);
    out.push('\n');
    for (field, value) in fields {
        out.push_str(field.as_ref());
        out.push_str(": ");
        out.push_str(value);
        out.push('\n');
    }
    out.push('\n');
}

fn format_vector_f64(values: &[f64]) -> String {
    let inner = values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{inner}]")
}

fn aot_toolchain_label(toolchain: Lsode2AotToolchain) -> &'static str {
    match toolchain {
        Lsode2AotToolchain::CTcc => "c_tcc",
        Lsode2AotToolchain::CGcc => "c_gcc",
        Lsode2AotToolchain::Zig => "zig",
        Lsode2AotToolchain::Rust => "rust",
    }
}

fn aot_profile_label(profile: Lsode2AotProfile) -> &'static str {
    match profile {
        Lsode2AotProfile::Debug => "debug",
        Lsode2AotProfile::Release => "release",
    }
}

fn linear_solver_choice_label(choice: Lsode2LinearSolverChoice) -> &'static str {
    match choice {
        Lsode2LinearSolverChoice::DenseLu => "dense_lu",
        Lsode2LinearSolverChoice::FaerSparseLu => "faer_sparse_lu",
        Lsode2LinearSolverChoice::LapackFaithfulBandedLu => "lapack_faithful_banded_lu",
    }
}

#[cfg(test)]
#[path = "task_parser_reactor_ivp_tests.rs"]
mod task_parser_reactor_ivp_tests;
