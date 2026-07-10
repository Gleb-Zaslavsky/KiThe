//! # Task Parser for Reactor BVP Module
//!
//! This module provides comprehensive parsing capabilities for chemical reactor boundary value problems (BVPs)
//! from structured text files. It bridges the gap between human-readable configuration files and the complex
//! data structures required by numerical BVP solvers.
//!
//! ## Purpose
//!
//! The module addresses the challenge of configuring complex reactor simulations through text files rather
//! than hardcoded parameters. It enables users to define complete reactor problems including:
//! - Chemical kinetics (reactions, Arrhenius parameters)
//! - Thermodynamic properties (heat capacity, thermal conductivity)
//! - Transport properties (diffusion coefficients)
//! - Boundary conditions and initial conditions
//! - Solver settings (tolerances, bounds, numerical methods)
//!
//! ## Boundary between KiThe and RustedSciThe
//!
//! ### Stays in KiThe
//! - `process_conditions` and other physical reactor inputs
//! - reactions, species, thermophysical data, and transport data
//! - boundary-condition structure that belongs to the reactor model
//! - `initial_guess`, because this is a reactor-BVP-specific initial-condition strategy
//!   rather than a generic numerical-solver setting
//!
//! ### Moves to RustedSciThe
//! - tolerance and bound parsing for the solver engine
//! - matrix backend selection (`Sparse`, `Banded`)
//! - symbolic backend selection (`ExprLegacy`, `AtomView`)
//! - AOT / lambdify execution policy and compiler selection
//! - solver iteration, logging, chunking, and other backend execution knobs
//! - solver-side task parsing that RustedSciThe now understands natively
//!
//! The goal is for KiThe to parse the physics and initial-condition strategy, while
//! RustedSciThe owns the numerical solver contract and all backend-specific options.
//!
//! ## Main Methods
//!
//! - **`set_reactor_params_from_hashmap()`**: Core parsing method that extracts reactor parameters from DocumentMap
//! - **`parse_tolerance_and_bounds()`**: Specialized parser for solver tolerance and bounds configurations
//! - **`parse_file()`**: File I/O wrapper that creates DocumentParser from file path
//! - **`parse_parameters_with_exact_names()`**: High-level parsing orchestrator
//! - **`solve_from_map()`**: Complete workflow from parsed data to solved BVP
//! - **`solve_from_file()`**: One-shot method: file → parsing → solving
//!
//!
//! ## Practical split rule
//!
//! If a setting changes the physics of the reactor task, it stays here.
//! If a setting changes how the NRBVP backend solves the task, it should move to RustedSciThe.
//!
//! ## Non-Obvious Code Features & Tips
//!
//! ### Variable Shadowing Prevention
//! The code carefully avoids variable name conflicts when parsing nested structures.
//! Example: `problem_name_val` instead of reusing `process_conditions` in if-let statements.
//!
//! ### DocumentMap Navigation Pattern
//! Uses consistent `.get().expect().clone().unwrap()[0].as_type()` pattern for extracting
//! values from the nested HashMap structure returned by the document parser.
//!
//! ### Groups Parsing Logic
//! The groups parsing implements a two-stage process: first check if groups are enabled,
//! then iterate through substances to build the nested HashMap structure.
//!
//! ### Vector Extraction Trick
//! For arrays like thermal_effects and reaction parameters, uses `.as_vector().unwrap().clone()`
//! to extract Vec<f64> directly from the parser's Value enum.
//!
//! ### Bounds Tuple Construction
//! Bounds parsing extracts 2-element vectors and converts them to tuples using array indexing:
//! `(bounds_vec[0], bounds_vec[1])` for efficient tuple creation.
//!
//! ### Error Handling Strategy
//! Uses `.expect()` for required parameters and `if let Some()` for optional ones,
//! providing clear error messages for missing required configuration.

use super::SimpleReactorBVP::{FastElemReact, SimpleReactorTask};
use crate::ReactorsBVP::reactor_BVP_utils::InitialConfig;
use crate::ReactorsBVP::reactor_BVP_utils::{BoundsConfig, ScalingConfig, ToleranceConfig};
use crate::Utils::show_this_pic::show_image;
use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
use log::info;
use nalgebra::DMatrix;
use std::collections::HashMap;
use RustedSciThe::numerical::BVP_Damp::task_parser_damped;

fn missing_field(section: &str, field: &str) -> crate::ReactorsBVP::SimpleReactorBVP::ReactorError {
    crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(format!(
        "Missing `{}` in `{}` section",
        field, section
    ))
}

fn invalid_numeric(
    section: &str,
    field: &str,
    value: f64,
) -> crate::ReactorsBVP::SimpleReactorBVP::ReactorError {
    crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidNumericValue(format!(
        "`{}.{} = {}` is not finite",
        section, field, value
    ))
}

fn required_section<'a>(
    task_hashmap: &'a DocumentMap,
    section: &str,
) -> Result<
    &'a HashMap<String, Option<Vec<Value>>>,
    crate::ReactorsBVP::SimpleReactorBVP::ReactorError,
> {
    task_hashmap
        .get(section)
        .ok_or_else(|| missing_field("document", section))
}

fn required_values<'a>(
    section: &'a HashMap<String, Option<Vec<Value>>>,
    section_name: &str,
    field: &str,
) -> Result<&'a [Value], crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    section
        .get(field)
        .ok_or_else(|| missing_field(section_name, field))?
        .as_ref()
        .map(|values| values.as_slice())
        .ok_or_else(|| missing_field(section_name, field))
}

fn required_f64(
    section: &HashMap<String, Option<Vec<Value>>>,
    section_name: &str,
    field: &str,
) -> Result<f64, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let value = required_values(section, section_name, field)?
        .first()
        .ok_or_else(|| missing_field(section_name, field))?
        .as_float()
        .ok_or_else(|| missing_field(section_name, field))?;
    if !value.is_finite() {
        return Err(invalid_numeric(section_name, field, value));
    }
    Ok(value)
}

fn required_usize(
    section: &HashMap<String, Option<Vec<Value>>>,
    section_name: &str,
    field: &str,
) -> Result<usize, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let value = required_values(section, section_name, field)?
        .first()
        .ok_or_else(|| missing_field(section_name, field))?;
    value_as_usize(value).ok_or_else(|| missing_field(section_name, field))
}

fn optional_string(
    section: &HashMap<String, Option<Vec<Value>>>,
    field: &str,
) -> Result<Option<String>, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    Ok(match section.get(field) {
        None => None,
        Some(values) => {
            let values = values
                .as_ref()
                .ok_or_else(|| missing_field("section", field))?;
            let value = values
                .first()
                .ok_or_else(|| missing_field("section", field))?;
            value
                .as_option_string()
                .cloned()
                .or_else(|| value.as_string().cloned())
        }
    })
}

// RustedSciThe's parser historically emits both Integer and Usize values.
// Keeping this conversion local avoids relying on one exact parser variant.
fn value_as_usize(value: &Value) -> Option<usize> {
    match value {
        Value::Usize(value) => Some(*value),
        Value::Integer(value) if *value >= 0 => Some(*value as usize),
        Value::Float(value) if value.is_finite() && *value >= 0.0 => Some(*value as usize),
        Value::String(value) => value.trim().parse::<usize>().ok(),
        Value::Optional(Some(inner)) => value_as_usize(inner),
        _ => None,
    }
}

fn required_vector_f64(
    section: &HashMap<String, Option<Vec<Value>>>,
    section_name: &str,
    field: &str,
) -> Result<Vec<f64>, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let values = required_values(section, section_name, field)?;
    let vector = values
        .first()
        .ok_or_else(|| missing_field(section_name, field))?
        .as_vector()
        .ok_or_else(|| missing_field(section_name, field))?;
    if let Some(invalid) = vector.iter().copied().find(|value| !value.is_finite()) {
        return Err(invalid_numeric(section_name, field, invalid));
    }
    Ok(vector.clone())
}

fn required_bool(
    section: &HashMap<String, Option<Vec<Value>>>,
    section_name: &str,
    field: &str,
) -> Result<bool, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    required_values(section, section_name, field)?
        .first()
        .ok_or_else(|| missing_field(section_name, field))?
        .as_boolean()
        .ok_or_else(|| missing_field(section_name, field))
}

// ============================= KiThe physics / run-setup block =============================
//
// This is still part of the reactor task contract. It defines the physical run window,
// not the numerical method used to solve it.
fn parse_physics_run_settings_from_map(
    task_hashmap: &DocumentMap,
) -> Result<(f64, f64, usize, String), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let process_conditions = required_section(task_hashmap, "process_conditions")?;
    let t0 = required_f64(process_conditions, "process_conditions", "t0")?;
    let t_end = required_f64(process_conditions, "process_conditions", "t_end")?;
    let n_steps = required_usize(process_conditions, "process_conditions", "n_steps")?;
    let arg = match process_conditions.get("arg") {
        Some(value) => {
            let values = value
                .as_ref()
                .ok_or_else(|| missing_field("process_conditions", "arg"))?;
            let arg = values
                .first()
                .ok_or_else(|| missing_field("process_conditions", "arg"))?
                .as_string()
                .cloned()
                .ok_or_else(|| missing_field("process_conditions", "arg"))?;
            Some(arg)
        }
        None => None,
    }
    .unwrap_or_else(|| "x".to_string());
    Ok((t0, t_end, n_steps, arg))
}

fn set_initial_guess_from_map_data(
    task_hashmap: &DocumentMap,
    n_steps: usize,
    unknown_count: usize,
) -> Result<DMatrix<f64>, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let inconfing = InitialConfig::new();
    let initial_guess: DMatrix<f64> = if let Some(guess_info) = task_hashmap.get("initial_guess") {
        let res = if let Some(universal) = guess_info.get("universal") {
            let values = universal
                .as_ref()
                .ok_or_else(|| missing_field("initial_guess", "universal"))?;
            let val = values
                .first()
                .ok_or_else(|| missing_field("initial_guess", "universal"))?
                .as_float()
                .ok_or_else(|| missing_field("initial_guess", "universal"))?;
            if !val.is_finite() {
                return Err(invalid_numeric("initial_guess", "universal", val));
            }
            inconfing.only_one_value_for_all_initial(val, n_steps, unknown_count)?
        } else {
            let C_val = required_f64(guess_info, "initial_guess", "C")?;
            let J_val = required_f64(guess_info, "initial_guess", "J")?;
            let Teta_val = required_f64(guess_info, "initial_guess", "Teta")?;
            let q_val = required_f64(guess_info, "initial_guess", "q")?;
            let map = HashMap::from([
                ("C".to_string(), C_val),
                ("J".to_string(), J_val),
                ("Teta".to_string(), Teta_val),
                ("q".to_string(), q_val),
            ]);
            inconfing.all_const_initial(map, n_steps)?
        };
        res
    } else {
        inconfing.only_one_value_for_all_initial(1e-2, n_steps, unknown_count)?
    };
    Ok(initial_guess)
}

/// Normalize solver-facing compatibility keys before delegating to RustedSciThe.
///
/// KiThe keeps the user-facing names readable, while the native damped parser
/// expects a small set of typed fields and canonical spellings.
pub(crate) fn normalize_reactor_physics_task_map(task_map: &DocumentMap) -> DocumentMap {
    let mut normalized = task_map.clone();
    normalize_exclusive_solver_backend(&mut normalized);
    if let Some(grid_refinement) = normalized.get_mut("grid_refinement") {
        if let Some(values) = grid_refinement.remove("grcar_smooke") {
            grid_refinement
                .entry("grcarsmooke".to_string())
                .or_insert(values);
        }
        if let Some(values) = grid_refinement.remove("double_points") {
            grid_refinement
                .entry("doubleoints".to_string())
                .or_insert(values);
        }
    }
    if adaptive_grid_is_explicitly_disabled(&normalized) {
        // RustedSciThe treats the presence of `adaptive_strategy` as an enabled
        // adaptive grid contract. A user-level `strategy_params.adaptive: None`
        // is the opposite: keep damping settings, but disable all grid adaptation.
        normalized.remove("adaptive_strategy");
        normalized.remove("grid_refinement");
    }
    normalize_solver_counter_field(&mut normalized, "solver_settings", "max_iterations");
    normalize_solver_counter_field(&mut normalized, "solver_settings", "refinement_steps");
    normalize_solver_counter_field(&mut normalized, "strategy_params", "max_jac");
    normalize_solver_counter_field(&mut normalized, "strategy_params", "max_damp_iter");
    normalize_solver_counter_field(&mut normalized, "adaptive_strategy", "version");
    normalize_solver_counter_field(&mut normalized, "adaptive_strategy", "max_refinements");
    normalized
}

const AOT_ONLY_SOLVER_FIELDS: &[&str] = &[
    "aot_codegen_backend",
    "aot_c_compiler",
    "aot_build_policy",
    "aot_build_profile",
    "aot_compile_preset",
    "aot_execution_policy",
];

#[derive(Clone, Copy)]
enum ExclusiveSolverBackend {
    Lambdify,
    Aot,
}

/// Normalize KiThe's exclusive execution choice before invoking RustedSciThe.
///
/// RustedSciThe permits artifact lifecycle and callback selection to be
/// configured independently. KiThe intentionally exposes Lambdify/AOT as an
/// exclusive choice, so old documents containing `lambdify_only` together with
/// `build_if_missing` are reduced to one coherent backend contract here.
fn normalize_exclusive_solver_backend(document: &mut DocumentMap) {
    let Some(section) = document.get_mut("solver_settings") else {
        return;
    };

    let generated_backend = section
        .get("generated_backend")
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(Value::as_string)
        .map(|value| value.trim().to_ascii_lowercase());
    let backend_policy = section
        .get("backend_policy")
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(Value::as_string)
        .map(|value| value.trim().to_ascii_lowercase());

    let generated_backend_mode = generated_backend
        .as_deref()
        .and_then(|backend| {
            if backend.contains("_aot") || backend.starts_with("aot_") {
                Some(ExclusiveSolverBackend::Aot)
            } else if backend.contains("lambdify") {
                Some(ExclusiveSolverBackend::Lambdify)
            } else {
                None
            }
        });
    // Preserve explicit unknown presets so the native typed parser can reject
    // them instead of silently changing user intent.
    if generated_backend.is_some() && generated_backend_mode.is_none() {
        return;
    }

    let backend = generated_backend_mode
        .or_else(|| match backend_policy.as_deref() {
            Some("aot_only") => Some(ExclusiveSolverBackend::Aot),
            Some("lambdify_only") => Some(ExclusiveSolverBackend::Lambdify),
            _ => None,
        })
        .unwrap_or(ExclusiveSolverBackend::Lambdify);

    let method = section
        .get("method")
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(Value::as_string)
        .map(|value| value.trim().to_ascii_lowercase())
        .unwrap_or_else(|| "banded".to_string());
    let matrix_prefix = if method == "sparse" { "sparse" } else { "banded" };

    match backend {
        ExclusiveSolverBackend::Lambdify => {
            section.insert(
                "generated_backend".to_string(),
                Some(vec![Value::String(format!("{matrix_prefix}_lambdify"))]),
            );
            section.insert(
                "backend_policy".to_string(),
                Some(vec![Value::String("lambdify_only".to_string())]),
            );
            for field_name in AOT_ONLY_SOLVER_FIELDS {
                section.remove(*field_name);
            }
        }
        ExclusiveSolverBackend::Aot => {
            if generated_backend.is_none() {
                let compiler = section
                    .get("aot_c_compiler")
                    .and_then(|slot| slot.as_ref())
                    .and_then(|values| values.first())
                    .and_then(Value::as_string)
                    .map(|value| value.trim().to_ascii_lowercase())
                    .unwrap_or_else(|| "tcc".to_string());
                section.insert(
                    "generated_backend".to_string(),
                    Some(vec![Value::String(format!(
                        "{matrix_prefix}_aot_{compiler}"
                    ))]),
                );
            }
            section.insert(
                "backend_policy".to_string(),
                Some(vec![Value::String("aot_only".to_string())]),
            );
        }
    }
}

/// Returns true when the task explicitly disables adaptive grid refinement.
///
/// The old and still-valid reactor task syntax keeps this switch in
/// `strategy_params` as `adaptive: None`. GUI rendering can create empty
/// `adaptive_strategy` / `grid_refinement` sections, so normalization must let
/// the explicit switch win before delegating to the native RustedSciThe parser.
fn adaptive_grid_is_explicitly_disabled(document: &DocumentMap) -> bool {
    document
        .get("strategy_params")
        .and_then(|section| section.get("adaptive"))
        .is_some_and(|slot| slot_is_explicit_none(slot.as_ref()))
}

/// Checks whether a document field carries an explicit `None` marker.
fn slot_is_explicit_none(values: Option<&Vec<Value>>) -> bool {
    values
        .and_then(|values| values.first())
        .is_some_and(value_is_explicit_none)
}

/// Checks the value shapes that can represent `None` after GUI edits or parsing.
fn value_is_explicit_none(value: &Value) -> bool {
    match value {
        Value::Optional(None) => true,
        Value::String(value) => value.trim().eq_ignore_ascii_case("none"),
        _ => false,
    }
}

/// Coerce a solver counter field into the scalar shape accepted by RustedSciThe.
///
/// RustedSciThe 0.4.12 names these fields as `usize`, but its `Value::as_usize`
/// helper currently reads non-negative `Value::Integer` payloads. KiThe keeps
/// the GUI flexible and normalizes to that native parser contract at the handoff.
fn normalize_solver_counter_field(
    document: &mut DocumentMap,
    section_name: &str,
    field_name: &str,
) {
    let Some(section) = document.get_mut(section_name) else {
        return;
    };
    let Some(values) = section.get_mut(field_name).and_then(|slot| slot.as_mut()) else {
        return;
    };

    if let Some(first) = values.first_mut() {
        if let Some(value) = value_as_usize(first) {
            *first = Value::Integer(value as i64);
            values.truncate(1);
        }
    }
}

/// Convert the typed damped solver-family settings into the canonical tolerance map.
///
/// The native RustedSciThe parser owns the typed solver contract. KiThe still needs
/// the expanded `C/J/Teta/q` family because the reactor solver stores tolerances in
/// the per-unknown map format used by the physical model.
fn tolerance_config_from_damped_spec(
    spec: &task_parser_damped::BvpDampedSolverSettingsSpec,
) -> Result<ToleranceConfig, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let rel_tolerance = spec.rel_tolerance.as_ref().ok_or_else(|| {
        crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
            "Missing `rel_tolerance` in `solver_settings` section".to_string(),
        )
    })?;
    Ok(ToleranceConfig::new(
        *rel_tolerance
            .get("C")
            .ok_or_else(|| missing_field("rel_tolerance", "C"))?,
        *rel_tolerance
            .get("J")
            .ok_or_else(|| missing_field("rel_tolerance", "J"))?,
        *rel_tolerance
            .get("Teta")
            .ok_or_else(|| missing_field("rel_tolerance", "Teta"))?,
        *rel_tolerance
            .get("q")
            .ok_or_else(|| missing_field("rel_tolerance", "q"))?,
    ))
}

/// Convert the typed damped solver-family settings into the canonical bounds map.
///
/// This mirrors [`tolerance_config_from_damped_spec`] and keeps the boundary between
/// solver-document parsing and reactor-state handoff explicit.
fn bounds_config_from_damped_spec(
    spec: &task_parser_damped::BvpDampedSolverSettingsSpec,
) -> Result<BoundsConfig, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let bounds = spec.bounds.as_ref().ok_or_else(|| {
        crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
            "Missing `bounds` in `solver_settings` section".to_string(),
        )
    })?;
    Ok(BoundsConfig::new(
        *bounds
            .get("C")
            .ok_or_else(|| missing_field("bounds", "C"))?,
        *bounds
            .get("J")
            .ok_or_else(|| missing_field("bounds", "J"))?,
        *bounds
            .get("Teta")
            .ok_or_else(|| missing_field("bounds", "Teta"))?,
        *bounds
            .get("q")
            .ok_or_else(|| missing_field("bounds", "q"))?,
    ))
}

/// Validate the dense solution matrix produced by the solver.
///
/// Solver backends may report a result without surfacing a structured failure.
/// This helper keeps the public contract strict by rejecting empty or non-finite
/// matrices before they are stored in reactor state.
fn validate_solver_output_matrix(
    solution: &DMatrix<f64>,
) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    if solution.nrows() == 0 || solution.ncols() == 0 {
        return Err(
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::CalculationError(
                "Solver returned an empty solution matrix".to_string(),
            ),
        );
    }
    if let Some((idx, value)) = solution
        .iter()
        .copied()
        .enumerate()
        .find(|(_, value)| !value.is_finite())
    {
        return Err(
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidNumericValue(format!(
                "Solver returned non-finite value {} at flat index {}",
                value, idx
            )),
        );
    }
    Ok(())
}

/// Build, configure, solve, and store a validated NRBVP snapshot.
///
/// The parser-facing solve methods share the same handoff logic; keeping it in
/// one place avoids drift between `solve_from_map` and `solve_from_parsed`.
fn solve_and_store_nrbvp(
    reactor: &mut SimpleReactorTask,
    task_map: &DocumentMap,
    t0: f64,
    t_end: f64,
    n_steps: usize,
    arg: String,
    initial_guess: DMatrix<f64>,
) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    let solver_task_map = normalize_reactor_physics_task_map(task_map);
    let damped_spec = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&solver_task_map)
        .map_err(|err| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidConfiguration(
                format!("RustedSciThe rejected solver settings: {err}"),
            )
        })?;
    // The typed RST options already contain tolerances, bounds, and backend
    // execution details, so the KiThe handoff only needs to pass them through.
    let solver_options = damped_spec
        .build_solver_options()
        .map_err(|err| crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidConfiguration(
            format!("RustedSciThe rejected solver settings: {err}")
        ))?;
    // The backend still expects the reactor-specific expanded maps, so we keep
    // this small compatibility bridge until RustedSciThe accepts them directly.
    let tolerance_config = tolerance_config_from_damped_spec(&damped_spec)?;
    let bounds_config = bounds_config_from_damped_spec(&damped_spec)?;
    let full_bounds = Some(bounds_config.to_full_bounds_map(&reactor.kindata.substances));
    let full_rel_tolerance = Some(tolerance_config.to_full_tolerance_map(&reactor.kindata.substances));

    let solver = &reactor.solver;
    let mut nr = solver.build_nrbvp_backend(
        crate::ReactorsBVP::SimpleReactorBVP::NrbvpHandoffConfig::new(
            initial_guess,
            t0,
            t_end,
            n_steps,
            damped_spec.scheme.clone(),
            damped_spec.strategy.clone(),
            damped_spec.strategy_params.clone(),
            damped_spec.linear_sys_method.clone(),
            damped_spec.method.clone(),
            damped_spec.abs_tolerance,
            full_rel_tolerance,
            damped_spec.max_iterations,
            full_bounds,
            damped_spec.loglevel.clone(),
            damped_spec.dont_save_log,
        )
        .with_solver_options(solver_options),
    )?;
    nr.before_solve_preprocessing();
    nr.solve();

    let solution = nr.get_result().ok_or_else(|| {
        crate::ReactorsBVP::SimpleReactorBVP::ReactorError::CalculationError(
            "Solver finished without producing a solution matrix".to_string(),
        )
    })?;
    validate_solver_output_matrix(&solution)?;
    if nr.x_mesh.is_empty() {
        return Err(
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::CalculationError(
                "Solver finished without producing an x mesh".to_string(),
            ),
        );
    }
    if nr.x_mesh.iter().any(|value| !value.is_finite()) {
        return Err(
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidNumericValue(
                "Solver x mesh contains non-finite values".to_string(),
            ),
        );
    }

    reactor.solver.x_mesh = Some(nr.x_mesh);
    reactor.solver.solution = Some(solution);
    reactor.solver.arg_name = arg;
    set_postprocessing_from_map(reactor, task_map)?;
    Ok(())
}

/// Apply post-processing options from a parsed document map.
///
/// This keeps the post-processing stage independent from `DocumentParser`
/// so the solver path can work from either parsed or already-materialized
/// configuration data without reparsing or mutating parser state.
fn set_postprocessing_from_map(
    reactor: &mut SimpleReactorTask,
    task_hashmap: &DocumentMap,
) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    reactor.check_balances()?;
    let solver_settings = required_section(task_hashmap, "postprocessing")?;

    // Flag to create plot via Rust native crate.
    let plot_flag = match solver_settings.get("plot") {
        Some(_) => required_bool(solver_settings, "postprocessing", "plot")?,
        None => false,
    };
    // Flag to create plot via GNU plot library.
    let gnuplot_flag = match solver_settings.get("gnuplot") {
        Some(_) => required_bool(solver_settings, "postprocessing", "gnuplot")?,
        None => false,
    };
    // Flag to save solution to txt.
    let save_flag = match solver_settings.get("save") {
        Some(_) => required_bool(solver_settings, "postprocessing", "save")?,
        None => false,
    };
    let save_to_csv = match solver_settings.get("save_to_csv") {
        Some(_) => required_bool(solver_settings, "postprocessing", "save_to_csv")?,
        None => false,
    };

    let name = optional_string(solver_settings, "filename")?;
    // Return from dimensionless to dimensioned unknowns.
    let return_to_dimension = match solver_settings.get("return_to_dimension") {
        Some(_) => required_bool(solver_settings, "postprocessing", "return_to_dimension")?,
        None => true,
    };
    let no_plots_in_terminal = match solver_settings.get("no_plots_in_terminal") {
        Some(_) => required_bool(solver_settings, "postprocessing", "no_plots_in_terminal")?,
        None => false,
    };

    if return_to_dimension {
        reactor.postprocessing()?;
    }
    if plot_flag {
        reactor.plot()?;
        let _ = show_image("Teta");
    }
    if gnuplot_flag {
        reactor.gnuplot()?;
        let _ = show_image("Teta");
    }
    if !no_plots_in_terminal {
        reactor.plot_in_terminal()?;
    }
    if save_flag {
        reactor.save_to_file(name.clone())?
    }
    if save_to_csv {
        reactor.save_to_csv(name)?;
    }
    reactor.estimate_values()?;
    Ok(())
}
impl SimpleReactorTask {
    /// Parses reactor parameters from a DocumentMap and populates the SimpleReactorTask structure.
    ///
    /// This is the core parsing method that extracts all reactor configuration from the parsed
    /// document structure. It handles reactions, process conditions, boundary conditions,
    /// diffusion coefficients, and chemical groups.
    ///
    /// # Arguments
    /// * `task_hashmap` - The parsed document structure containing all configuration sections
    ///
    /// # Sections Parsed
    /// - `reactions`: Arrhenius parameters (A, n, E, Q) for each reaction
    /// - `process_conditions`: Physical properties (Tm, P, Cp, Lambda, m, etc.)
    /// - `boundary_condition`: Initial/boundary values for all variables
    /// - `diffusion_coefficients`: Transport properties for each substance
    /// - Chemical groups: Atomic composition for custom substance definitions
    pub fn set_reactor_params_from_hashmap(
        &mut self,
        task_hashmap: &DocumentMap,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let reactions = required_section(task_hashmap, "reactions")?;
        let mut vec_of_struct: Vec<FastElemReact> = Vec::with_capacity(reactions.len());
        for (key, value) in reactions.iter() {
            let values = value
                .as_ref()
                .ok_or_else(|| missing_field("reactions", key))?;
            let rates = values
                .first()
                .ok_or_else(|| missing_field("reactions", key))?
                .as_vector()
                .ok_or_else(|| missing_field("reactions", key))?;
            if rates.len() != 4 {
                return Err(
                    crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidConfiguration(
                        format!(
                            "Reaction `{}` must store exactly 4 Arrhenius parameters",
                            key
                        ),
                    ),
                );
            }
            if let Some(invalid) = rates.iter().copied().find(|value| !value.is_finite()) {
                return Err(invalid_numeric("reactions", key, invalid));
            }
            vec_of_struct.push(FastElemReact {
                eq: key.to_owned(),
                A: rates[0],
                n: rates[1],
                E: rates[2],
                Q: rates[3],
            });
        }
        self.fast_react_set(vec_of_struct)?;

        let process_conditions = required_section(task_hashmap, "process_conditions")?;
        self.problem_name = optional_string(process_conditions, "problem_name")?;
        self.problem_description = optional_string(process_conditions, "problem_description")?;

        let substances_values =
            required_values(process_conditions, "process_conditions", "substances")?;
        let mut substances: Vec<String> = Vec::with_capacity(substances_values.len());
        for val in substances_values.iter() {
            let substance = val
                .as_string()
                .cloned()
                .ok_or_else(|| missing_field("process_conditions", "substances"))?;
            substances.push(substance);
        }
        self.kindata.substances = substances;

        self.Tm = required_f64(process_conditions, "process_conditions", "Tm")?;
        let L = required_f64(process_conditions, "process_conditions", "L")?;
        let dT = required_f64(process_conditions, "process_conditions", "dT")?;
        let T_scale = required_f64(process_conditions, "process_conditions", "T_scale")?;
        self.scaling = ScalingConfig::new(dT, L, T_scale);
        self.P = required_f64(process_conditions, "process_conditions", "P")?;
        self.Cp = required_f64(process_conditions, "process_conditions", "Cp")?;
        self.Lambda = required_f64(process_conditions, "process_conditions", "Lambda")?;
        self.m = required_f64(process_conditions, "process_conditions", "m")?;
        if let Some(mass) = process_conditions.get("M") {
            let values = mass
                .as_ref()
                .ok_or_else(|| missing_field("process_conditions", "M"))?;
            let value = values
                .first()
                .ok_or_else(|| missing_field("process_conditions", "M"))?
                .as_float()
                .ok_or_else(|| missing_field("process_conditions", "M"))?;
            if !value.is_finite() {
                return Err(invalid_numeric("process_conditions", "M", value));
            }
            self.M = value;
        }

        self.thermal_effects =
            required_vector_f64(process_conditions, "process_conditions", "thermal_effects")?;

        let boundary_condition_section = required_section(task_hashmap, "boundary_condition")?;
        let mut boundary_condition: HashMap<String, f64> =
            HashMap::with_capacity(boundary_condition_section.len());
        for (key, value) in boundary_condition_section.iter() {
            let values = value
                .as_ref()
                .ok_or_else(|| missing_field("boundary_condition", key))?;
            let value = values
                .first()
                .ok_or_else(|| missing_field("boundary_condition", key))?
                .as_float()
                .ok_or_else(|| missing_field("boundary_condition", key))?;
            if !value.is_finite() {
                return Err(invalid_numeric("boundary_condition", key, value));
            }
            boundary_condition.insert(key.to_owned(), value);
        }
        self.boundary_condition = boundary_condition;

        let diffusion_coefficients_section =
            required_section(task_hashmap, "diffusion_coefficients")?;
        let mut diffusion_coefficients: HashMap<String, f64> =
            HashMap::with_capacity(diffusion_coefficients_section.len());
        for (key, value) in diffusion_coefficients_section.iter() {
            let values = value
                .as_ref()
                .ok_or_else(|| missing_field("diffusion_coefficients", key))?;
            let value = values
                .first()
                .ok_or_else(|| missing_field("diffusion_coefficients", key))?
                .as_float()
                .ok_or_else(|| missing_field("diffusion_coefficients", key))?;
            if !value.is_finite() {
                return Err(invalid_numeric("diffusion_coefficients", key, value));
            }
            diffusion_coefficients.insert(key.to_owned(), value);
        }
        self.Diffusion = diffusion_coefficients;

        if let Some(groups_val) = process_conditions.get("groups") {
            let values = groups_val
                .as_ref()
                .ok_or_else(|| missing_field("process_conditions", "groups"))?;
            let groups_enabled = values
                .first()
                .ok_or_else(|| missing_field("process_conditions", "groups"))?
                .as_boolean()
                .ok_or_else(|| missing_field("process_conditions", "groups"))?;
            if groups_enabled {
                let mut groups: HashMap<String, HashMap<String, usize>> =
                    HashMap::with_capacity(self.kindata.substances.len());
                for sub_i in self.kindata.substances.iter() {
                    if let Some(sub_i_groups) = task_hashmap.get(sub_i.as_str()) {
                        let mut sub_i_groups_hashmap: HashMap<String, usize> =
                            HashMap::with_capacity(sub_i_groups.len());
                        for (key, value) in sub_i_groups.iter() {
                            let values = value.as_ref().ok_or_else(|| missing_field(sub_i, key))?;
                            let value = values
                                .first()
                                .ok_or_else(|| missing_field(sub_i, key))?
                                .as_usize()
                                .ok_or_else(|| missing_field(sub_i, key))?;
                            sub_i_groups_hashmap.insert(key.to_owned(), value);
                        }
                        groups.insert(sub_i.clone(), sub_i_groups_hashmap);
                    }
                }
                self.kindata.groups = Some(groups);
            }
        }

        Ok(())
    }
    /// Parses solver tolerance and bounds configurations from the document.
    ///
    /// This is solver-side configuration, not reactor-physics parsing.
    ///
    /// Extracts numerical solver settings including relative tolerances and variable bounds
    /// for the four main variable types: C (concentrations), J (fluxes), Teta (temperature), q (heat flux).
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the parsed configuration
    ///
    /// # Returns
    /// * `Ok((ToleranceConfig, BoundsConfig))` - Parsed tolerance and bounds configurations
    /// * `Err(String)` - Error message if required sections are missing
    ///
    /// # Expected Document Structure
    /// ```text
    /// rel_tolerance
    /// C: 1e-5
    /// J: 1e-5
    /// Teta: 1e-5
    /// q: 1e-5
    /// bounds
    /// C: -10.0, 10.0
    /// J: -1e20, 1e20
    /// Teta: -100.0, 100.0
    /// q: -1e20, 1e20
    /// ```
    fn parse_tolerance_and_bounds(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(ToleranceConfig, BoundsConfig), crate::ReactorsBVP::SimpleReactorBVP::ReactorError>
    {
        let task_hashmap = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "document must be parsed before tolerance extraction".to_string(),
            )
        })?;
        let normalized = normalize_reactor_physics_task_map(task_hashmap);
        let spec = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&normalized)
            .map_err(|err| {
                crate::ReactorsBVP::SimpleReactorBVP::ReactorError::InvalidConfiguration(
                    format!("RustedSciThe rejected solver settings: {err}"),
                )
            })?;
        let tolerance_config = tolerance_config_from_damped_spec(&spec)?;
        let bounds_config = bounds_config_from_damped_spec(&spec)?;
        Ok((tolerance_config, bounds_config))
    }

    /// Legacy spelling kept for compatibility with older call sites.
    fn parse_toleranse_and_bounds(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(ToleranceConfig, BoundsConfig), crate::ReactorsBVP::SimpleReactorBVP::ReactorError>
    {
        self.parse_tolerance_and_bounds(parser)
    }

    /// Extracts the physical run window from the parsed document.
    ///
    /// These are reactor-task inputs, not numerical solver knobs.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the configuration
    ///
    /// # Returns
    /// * `(t0, t_end, n_steps)` - Start time, end time, and number of grid steps
    pub fn parse_basic_settings(
        &mut self,
        parser: DocumentParser,
    ) -> Result<(f64, f64, usize, String), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let task_hashmap = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "document must be parsed before settings extraction".to_string(),
            )
        })?;
        parse_physics_run_settings_from_map(task_hashmap)
    }

    pub fn set_postprocessing_from_hashmap(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let result = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "document must be parsed before postprocessing".to_string(),
            )
        })?;
        set_postprocessing_from_map(self, result)
    }

    /// Legacy spelling kept for compatibility with older call sites.
    pub fn set_postpocessing_from_hashmap(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        self.set_postprocessing_from_hashmap(parser)
    }

    pub fn set_initial_guess_from_map(
        &mut self,
        parser: DocumentParser,
        n_steps: usize,
    ) -> Result<DMatrix<f64>, crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let result = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "document must be parsed before initial guess extraction".to_string(),
            )
        })?;
        set_initial_guess_from_map_data(result, n_steps, self.solver.unknowns.len())
    }
    /// Complete workflow: parses configuration, sets up BVP, configures solver, and solves the problem.
    ///
    /// This is a high-level method that orchestrates the entire process from parsed document
    /// to solved BVP. It handles all intermediate steps including BVP setup, solver configuration,
    /// tolerance/bounds parsing, and postprocessing.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the complete problem configuration
    ///
    /// # Process Flow
    /// 1. Parse reactor parameters
    /// 2. Setup BVP equations and boundary conditions
    /// 3. Configure NRBVP solver with parsed settings
    /// 4. Parse tolerances and bounds
    /// 5. Set initial guess and solve
    /// 6. Apply postprocessing (plotting, data export)
    pub fn solve_from_map(
        &mut self,
        mut parser: DocumentParser,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        self.parse_parameters_with_exact_names(&mut parser)?;
        let task_map = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "document must be parsed before solver setup".to_string(),
            )
        })?;
        let (t0, t_end, n_steps, arg) = parse_physics_run_settings_from_map(task_map)?;
        self.setup_bvp()?;
        if self.solver.unknowns.is_empty() {
            return Err(
                crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
                    "Solver unknowns were not initialized".to_string(),
                ),
            );
        }
        if self.solver.eq_system.is_empty() {
            return Err(
                crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
                    "Solver equation system was not initialized".to_string(),
                ),
            );
        }
        let initial_guess =
            set_initial_guess_from_map_data(task_map, n_steps, self.solver.unknowns.len())?;
        solve_and_store_nrbvp(
            self,
            task_map,
            t0,
            t_end,
            n_steps,
            arg,
            initial_guess,
        )
    }

    /// One-shot method: reads configuration file, parses it, and solves the BVP problem.
    ///
    /// Convenience method that combines file reading, parsing, and solving in a single call.
    /// Ideal for batch processing or when the entire workflow can be automated.
    ///
    /// # Arguments
    /// * `path` - Path to the configuration file
    pub fn solve_from_file(
        &mut self,
        path: std::path::PathBuf,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let mut parser = self
            .parse_file(Some(path))
            .map_err(crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError)?;
        parser
            .parse_document()
            .map_err(crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError)?;
        self.solve_from_map(parser)
    }

    /// Creates a DocumentParser from a file path.
    ///
    /// Simple wrapper around DocumentParser creation that handles file I/O.
    ///
    /// # Arguments
    /// * `path` - Optional path to configuration file
    ///
    /// # Returns
    /// * `Ok(DocumentParser)` - Configured parser ready for document parsing
    /// * `Err(String)` - File I/O error message
    pub fn parse_file(
        &mut self,
        path: Option<std::path::PathBuf>,
    ) -> Result<DocumentParser, String> {
        let mut parser = DocumentParser::new(String::new());
        parser.setting_from_file(path)?;
        Ok(parser)
    }

    /// High-level parsing orchestrator that processes the document and populates reactor parameters.
    ///
    /// Coordinates the document parsing process and delegates to `set_reactor_params_from_hashmap`
    /// for actual parameter extraction.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser with loaded configuration
    ///
    /// # Returns
    /// * `Ok(())` - Successful parsing
    /// * `Err(String)` - Parsing error message
    pub fn parse_parameters_with_exact_names(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        // Keep the helper idempotent: callers may hand us a parser that is already parsed.
        // Re-parsing can mutate the internal document shape in some backends, so we only
        // parse when the result cache is still empty.
        if parser.get_result().is_none() {
            parser
                .parse_document()
                .map_err(crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError)?;
        }
        let result = parser.get_result().ok_or_else(|| {
            crate::ReactorsBVP::SimpleReactorBVP::ReactorError::ParseError(
                "No result after parsing".to_string(),
            )
        })?;
        self.set_reactor_params_from_hashmap(result)?;
        Ok(())
    }

    pub fn solve_from_parsed(
        &mut self,
        map: DocumentMap,
    ) -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
        let task_map = map;
        self.set_reactor_params_from_hashmap(&task_map)?;
        let (t0, t_end, n_steps, arg) = parse_physics_run_settings_from_map(&task_map)?;
        self.setup_bvp()?;
        if self.solver.unknowns.is_empty() {
            return Err(
                crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
                    "Solver unknowns were not initialized".to_string(),
                ),
            );
        }
        if self.solver.eq_system.is_empty() {
            return Err(
                crate::ReactorsBVP::SimpleReactorBVP::ReactorError::MissingData(
                    "Solver equation system was not initialized".to_string(),
                ),
            );
        }
        let initial_guess =
            set_initial_guess_from_map_data(&task_map, n_steps, self.solver.unknowns.len())?;
        solve_and_store_nrbvp(
            self,
            &task_map,
            t0,
            t_end,
            n_steps,
            arg,
            initial_guess,
        )
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a template configuration file for reactor BVP tasks.
///
/// Generates a commented template file that users can fill with their specific values.
/// The template includes all required sections with example values and explanatory comments.
///
pub const SIMPLE_BVP_TEMPLATE: &'static str = r#"
        initial_guess
        universal:1e-2
        process_conditions
        problem_name: Some(some_name)
        problem_description: Some("some decription")
        substances: A, B
        t0: 0.0
        t_end: 1.0
        n_steps: 25
        arg:x
        Tm: 1500.0
        L: 5e-4
        dT: 600.0
        T_scale: 600.0
        P: 1e6
        Cp: 1464.4
        Lambda: 0.07
        m: 0.0043
        M: 0.0342
        thermal_effects: [102000.0]
        groups:true
        boundary_condition
        A: 0.999
        B: 0.001
        T: 800.0
        diffusion_coefficients
        A: 0.000009296
        B: 0.000009296
        group1
        H: 4
        N: 8
        C: 8
        O: 8
        group2
        H: 6
        C: 1
        O: 1
        reactions
        A=>10B: [130000.0, 0.0, 20920.0, 102000.0]
        solver_settings
        scheme: forward
        method: Banded
        strategy: Damped
        linear_sys_method: None
        generated_backend: banded_lambdify
        matrix_backend: Banded
        backend_policy: lambdify_only
        symbolic_backend: AtomView
        banded_linear_solver: auto
        refinement_steps: 5
        abs_tolerance: 1e-7
        max_iterations: 100
        loglevel: Some(info)
        dont_save_log: true
        bounds
        C: -10.0, 10.0
        J:  -1e20, 1e20
        Teta:-100.0, 100.0
        q: -1e20, 1e20
        rel_tolerance
        C: 1e-7
        J: 1e-7
        Teta: 1e-7
        q:  1e-7
        strategy_params
        max_jac: Some(3)
        max_damp_iter: Some(10)
        damp_factor: Some(0.5)
        # Adaptive grid refinement settings (optional)
        adaptive_strategy
        # Refinement version
        version: 1
        # Maximum refinement iterations
        max_refinements: 3
                
        #Grid refinement method and parameters
        grid_refinement
        // Available methods:
        // double_points: []
        // easy: [parameter]
        // grcar_smooke: [param1, param2, param3] (legacy token: grcarsmooke)
        // pearson: [param1, param2]
        // twopnt: [param1, param2, param3]
        grcar_smooke: [0.05, 0.05, 1.25]
        postprocessing
        gnuplot:true
        save_to_csv:false
        filename: meow
        gui_plot: true
        "#;
pub fn create_template() -> Result<(), crate::ReactorsBVP::SimpleReactorBVP::ReactorError> {
    use std::fs::File;
    use std::io::Write;

    let template_content = r#"
# Reactor BVP Configuration Template
// This template provides a starting point for configuring reactor BVP tasks.
// Users should fill in the values below for their specific problem.
// The template includes all required sections with example values and explanatory comments.
//
// some fields are filled with Some(value) it means the value in this field
// will be used in the task, otherwise the default value will be used
// if the field is None, the default value will be used
//
# Fill in the values below for your specific problem
# initial guess - if not set, default is 1e-2 for all variables and mesh
# if set universal:some_value - for all variables and mesh will be set this value
# if set C: some_value, J: some_value, Teta: some_value, q: some_value - for all 
# corresponing variables will be set this values
initial_guess:
universal: 1e-2
# Process conditions - main problem parameters
process_conditions
# Optional problem identification
problem_name: Some(YourProblemName)
problem_description: Some(YourProblemDescription)
# List of chemical substances (comma-separated)
substances: Substance1, Substance2
# Start position (real length is set in L field, so 
# t0 and t_end are dimensionless and reasonable to leave 
# them 0.0 and 1.0 )
t0: 0.0
# End position 
t_end: 1.0
# Number of grid points - the more you set 
# the more accurate the solution will be, but 
# the more time it will take to solve the problem
n_steps: 30
# Independent variable name (x, t, etc.)
arg: x
# Maximum temperature in the front of reaction [K]
Tm: 1500.0
# Characteristic length [m] - should be chosen in such 
# way that all reactions are ended at the end of the distance
L: 9e-4
# Temperature difference [K]
dT: 600.0
# Temperature scaling [K]
T_scale: 600.0
# Pressure [Pa]
P: 1e6
# Heat capacity [J/kg/K]
Cp: 1464.4
# Thermal conductivity [W/m/K]
Lambda: 0.07
# Mass velocity [kg/m2*s]
m: 0.0043
# Molar mass [kg/mol]
M: 0.0342
# Thermal effects for each reaction [J/mol]
thermal_effects: [102000.0]
# Enable custom atomic groups (true/false)
groups: true

# Boundary conditions - initial values for all variables
boundary_condition
# Initial concentration/mole fraction
Substance1: 0.999
# Initial concentration/mole fraction
Substance2: 0.001
# Initial temperature [K]
T: 800.0

# Diffusion coefficients for each substance [m²/s]
diffusion_coefficients
Substance1: 0.000009296
Substance2: 0.000009296

# Atomic composition (only if groups: true)
# Define atomic composition for each substance
Substance1
# Number of hydrogen atoms
H: 4
# Number of nitrogen atoms
N: 8
# Number of carbon atoms
C: 8
# Number of oxygen atoms
O: 8

Substance2
H: 6
C: 1
O: 1

# Chemical reactions with Arrhenius parameters
# Format: Reactant=>Product: [A, n, E, Q]
# A: pre-exponential factor, n: temperature exponent
# E: activation energy [J/mol/K],
# Q: heat of reaction [J/mol]
reactions
Substance1=>Substance2: [130000.0, 0.0, 20920.0, 102000.0]

# Numerical solver settings
solver_settings
# Discretization scheme
scheme: forward
# Solution method
method: Banded
# Convergence strategy
strategy: Damped
# Linear system solver
linear_sys_method: None
# High-level generated backend preset
generated_backend: banded_lambdify
# Generated matrix backend
matrix_backend: Banded
# Backend selection policy
backend_policy: lambdify_only
# Symbolic assembly backend
symbolic_backend: AtomView
# Banded linear solver policy
banded_linear_solver: auto
# Iterative refinement steps
refinement_steps: 5
# Absolute tolerance
abs_tolerance: 1e-5
# Maximum iterations
max_iterations: 100
# Logging level
loglevel: Some(info)
# Disable log saving
dont_save_log: true

# Variable bounds for solver
bounds
# Concentration bounds
C: -10.0, 10.0
# Flux bounds
J: -1e20, 1e20
# Temperature bounds
Teta: -100.0, 100.0
# Heat flux bounds
q: -1e20, 1e20

# Relative tolerances for each variable type
rel_tolerance
# Concentration tolerance
C: 1e-5
# Flux tolerance
J: 1e-5
# Temperature tolerance
Teta: 1e-5
# Heat flux tolerance
q: 1e-5

# Advanced solver strategy parameters
strategy_params
# Maximum Jacobian updates
max_jac: Some(3)
# Maximum damping iterations
max_damp_iter: Some(10)
# Damping factor
damp_factor: Some(0.5)

# Adaptive grid refinement settings (optional)
adaptive_strategy
# Refinement version
version: 1
# Maximum refinement iterations
max_refinements: 3
        
#Grid refinement method and parameters
grid_refinement
// Available methods:
// double_points: []
// easy: [parameter]
// grcar_smooke: [param1, param2, param3] (legacy token: grcarsmooke)
// pearson: [param1, param2]
// twopnt: [param1, param2, param3]
pearson: [0.1, 1.5]

# Output and visualization options
postprocessing
# Generate gnuplot output
gnuplot: true
# Save results to CSV
save_to_csv: false
# Output filename prefix
filename: output_name
"#;

    let mut file = File::create("template.txt").map_err(|err| {
        crate::ReactorsBVP::SimpleReactorBVP::ReactorError::CalculationError(format!(
            "Failed to create template.txt: {}",
            err
        ))
    })?;
    file.write_all(template_content.as_bytes()).map_err(|err| {
        crate::ReactorsBVP::SimpleReactorBVP::ReactorError::CalculationError(format!(
            "Failed to write template.txt: {}",
            err
        ))
    })?;
    info!("Template created: template.txt");
    Ok(())
}
////////////////////////////////////////////////////TESTS///////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;
    const task_content: &str = r#"
        initial_guess
        universal:1e-2
        process_conditions
        problem_name: Some(HMXTest)
        problem_description: Some(HMXdecompositiontest)
        substances: HMX, HMXprod
        t0: 0.0
        t_end: 1.0
        n_steps: 25
        arg:x
        Tm: 1500.0
        L: 5e-4
        dT: 600.0
        T_scale: 600.0
        P: 1e6
        Cp: 1464.4
        Lambda: 0.07
        m: 0.0043
        M: 0.0342
        thermal_effects: [102000.0]
        groups:true
        boundary_condition
        HMX: 0.999
        HMXprod: 0.001
        T: 800.0
        diffusion_coefficients
        HMX: 0.000009296
        HMXprod: 0.000009296
        HMX
        H: 4
        N: 8
        C: 8
        O: 8
        HMXprod
        H: 6
        C: 1
        O: 1
        reactions
        HMX=>10HMXprod: [130000.0, 0.0, 20920.0, 102000.0]
        solver_settings
        scheme: forward
        method: Banded
        strategy: Damped
        linear_sys_method: None
        generated_backend: banded_lambdify
        matrix_backend: Banded
        backend_policy: lambdify_only
        symbolic_backend: AtomView
        banded_linear_solver: auto
        refinement_steps: 5
        abs_tolerance: 1e-7
        max_iterations: 100
        loglevel: Some(info)
        dont_save_log: true
        bounds
        C: -10.0, 10.0
        J:  -1e20, 1e20
        Teta:-100.0, 100.0
        q: -1e20, 1e20
        rel_tolerance
        C: 1e-7
        J: 1e-7
        Teta: 1e-7
        q:  1e-7
        strategy_params
        max_jac: Some(3)
        max_damp_iter: Some(10)
        damp_factor: Some(0.5)
        # Adaptive grid refinement settings (optional)
        adaptive_strategy
        # Refinement version
        version: 1
        # Maximum refinement iterations
        max_refinements: 3
                
        #Grid refinement method and parameters
        grid_refinement
        // Available methods:
        // double_points: []
        // easy: [parameter]
        // grcar_smooke: [param1, param2, param3] (legacy token: grcarsmooke)
        // pearson: [param1, param2]
        // twopnt: [param1, param2, param3]
        grcar_smooke: [0.05, 0.05, 1.25]
        postprocessing
        gnuplot:true
        save_to_csv:false
        filename: meow
        "#;
    #[test]
    fn parse_directly() {
        let mut parser = DocumentParser::new(task_content.to_string());
        let result: DocumentMap = parser.parse_document().unwrap().clone();
        info!("result {:?}", result);
    }
    #[test]
    fn check_substances() {
        let mut parser = DocumentParser::new(task_content.to_string());
        let result: DocumentMap = parser.parse_document().unwrap().clone();
        let process_conditions = result.get("process_conditions").unwrap();
        let substances = process_conditions
            .get("substances")
            .unwrap()
            .clone()
            .unwrap();

        info!("substances {:?}", substances[0]);
    }

    #[test]
    fn test_solve_from_parsed() {
        let mut parser = DocumentParser::new(task_content.to_string());
        let result: DocumentMap = parser.parse_document().unwrap().clone();
        let mut reactor = SimpleReactorTask::new();
        reactor
            .solve_from_parsed(result)
            .expect("Failed to solve parsed task");
        assert!(reactor.solver.solution.is_some());
        assert!(reactor.solver.x_mesh.is_some());
        assert!(!reactor.solver.x_mesh.as_ref().unwrap().is_empty());
        info!("solution {:?}", reactor.solver.solution);
    }

    #[test]
    fn test_solve_from_parsed_with_explicit_postprocessing_flags() {
        let mut parser = DocumentParser::new(task_content.to_string());
        let mut result: DocumentMap = parser.parse_document().unwrap().clone();

        // Keep the regression path focused on the solve/postprocessing handoff
        // without creating file or terminal side effects during the test.
        let postprocessing = result
            .get_mut("postprocessing")
            .expect("postprocessing section must exist in the parsed task");
        postprocessing.insert(
            "plot".to_string(),
            Some(vec![Value::Boolean(false)]),
        );
        postprocessing.insert(
            "gnuplot".to_string(),
            Some(vec![Value::Boolean(false)]),
        );
        postprocessing.insert("save".to_string(), Some(vec![Value::Boolean(false)]));
        postprocessing.insert(
            "save_to_csv".to_string(),
            Some(vec![Value::Boolean(false)]),
        );
        postprocessing.insert(
            "return_to_dimension".to_string(),
            Some(vec![Value::Boolean(false)]),
        );
        postprocessing.insert(
            "no_plots_in_terminal".to_string(),
            Some(vec![Value::Boolean(true)]),
        );

        let mut reactor = SimpleReactorTask::new();
        reactor
            .solve_from_parsed(result)
            .expect("Failed to solve parsed task with explicit postprocessing flags");
        assert!(reactor.solver.solution.is_some());
        assert!(reactor.solver.x_mesh.is_some());
        assert!(!reactor.solver.x_mesh.as_ref().unwrap().is_empty());
        assert!(!reactor.solver.solution.as_ref().unwrap().is_empty());
    }
    #[test]
    fn parsng_task_elementary() {
        let mut reactor = SimpleReactorTask::new();

        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");

        // HMX data from create_hmx function

        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        info!("parser {:?}", parser);

        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("Failed to parse parameters");

        // Verify parsed data
        assert_eq!(reactor.problem_name, Some("HMXTest".to_string()));
        assert_eq!(
            reactor.kindata.substances,
            vec!["HMX".to_string(), "HMXprod".to_string()]
        );

        assert_eq!(reactor.Tm, 1500.0);
        assert_eq!(reactor.P, 1e6);
        assert_eq!(reactor.Cp, 1464.4);
        assert_eq!(reactor.Lambda, 0.07);
        assert_eq!(reactor.m, 0.0043);
        assert_eq!(reactor.M, 0.0342);
        assert_eq!(reactor.thermal_effects, vec![102000.0]);

        // Verify boundary conditions
        assert_eq!(reactor.boundary_condition.get("HMX"), Some(&0.999));
        assert_eq!(reactor.boundary_condition.get("HMXprod"), Some(&0.001));
        assert_eq!(reactor.boundary_condition.get("T"), Some(&800.0));

        // Verify diffusion coefficients
        assert_eq!(reactor.Diffusion.get("HMX"), Some(&0.000009296));
        assert_eq!(reactor.Diffusion.get("HMXprod"), Some(&0.000009296));

        // Verify groups
        assert!(reactor.kindata.groups.is_some());
        let groups = reactor.kindata.groups.as_ref().unwrap();
        assert_eq!(groups.get("HMX").unwrap().get("H"), Some(&4));
        assert_eq!(groups.get("HMX").unwrap().get("N"), Some(&8));
        assert_eq!(groups.get("HMXprod").unwrap().get("H"), Some(&6));
        assert_eq!(groups.get("HMXprod").unwrap().get("C"), Some(&1));

        // Verify reactions were parsed
        assert_eq!(reactor.kindata.vec_of_equations.len(), 1);
        assert_eq!(reactor.kindata.vec_of_equations[0], "HMX=>10HMXprod");
    }

    #[test]
    fn test_parse_toleranse_and_bounds() {
        let mut reactor = SimpleReactorTask::new();

        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test_bounds.txt");

        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Parse file
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();

        // Test parsing tolerances and bounds
        let (tolerance_config, bounds_config) = reactor
            .parse_toleranse_and_bounds(&mut parser)
            .expect("Failed to parse tolerances and bounds");

        // Verify tolerance config
        assert_eq!(tolerance_config.C, 1e-7);
        assert_eq!(tolerance_config.J, 1e-7);
        assert_eq!(tolerance_config.Teta, 1e-7);
        assert_eq!(tolerance_config.q, 1e-7);

        // Verify bounds config
        assert_eq!(bounds_config.C, (-10.0, 10.0));
        assert_eq!(bounds_config.J, (-1e20, 1e20));
        assert_eq!(bounds_config.Teta, (-100.0, 100.0));
        assert_eq!(bounds_config.q, (-1e20, 1e20));
    }

    #[test]
    fn test_build_bvp_damped_solver_options_accepts_aot_backend_options() {
        let mut solver_settings = HashMap::new();
        solver_settings.insert(
            "scheme".to_string(),
            Some(vec![Value::String("forward".to_string())]),
        );
        solver_settings.insert(
            "strategy".to_string(),
            Some(vec![Value::String("Damped".to_string())]),
        );
        solver_settings.insert(
            "method".to_string(),
            Some(vec![Value::String("Banded".to_string())]),
        );
        solver_settings.insert("abs_tolerance".to_string(), Some(vec![Value::Float(1e-8)]));
        solver_settings.insert("max_iterations".to_string(), Some(vec![Value::Integer(50)]));
        solver_settings.insert(
            "generated_backend".to_string(),
            Some(vec![Value::String("banded_aot_tcc".to_string())]),
        );
        solver_settings.insert(
            "matrix_backend".to_string(),
            Some(vec![Value::String("Banded".to_string())]),
        );
        solver_settings.insert(
            "symbolic_backend".to_string(),
            Some(vec![Value::String("AtomView".to_string())]),
        );
        solver_settings.insert(
            "aot_c_compiler".to_string(),
            Some(vec![Value::String("zig".to_string())]),
        );
        solver_settings.insert(
            "aot_build_policy".to_string(),
            Some(vec![Value::String("require_prebuilt".to_string())]),
        );
        solver_settings.insert(
            "aot_compile_preset".to_string(),
            Some(vec![Value::String("fast_build".to_string())]),
        );
        solver_settings.insert(
            "aot_execution_policy".to_string(),
            Some(vec![Value::String("auto".to_string())]),
        );

        let mut document = DocumentMap::new();
        document.insert("solver_settings".to_string(), solver_settings);

        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_settings = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&normalized)
            .expect("BVP damped solver settings should parse through the typed RST API");
        let options = solver_settings
            .build_solver_options()
            .expect("Typed damped settings should build RustedSciThe options");

        assert_eq!(
            options.generated_backend_config.aot_c_compiler.as_deref(),
            Some("zig")
        );
        assert!(matches!(
            options.generated_backend_config.aot_build_policy,
            RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::AotBuildPolicy::RequirePrebuilt
        ));
        assert!(matches!(
            options.generated_backend_config.aot_execution_policy,
            RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::AotExecutionPolicy::Auto
        ));
        assert!(matches!(
            options.generated_backend_config.backend_policy_override,
            Some(
                RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy::AotOnly
            )
        ));
        assert!(options.generated_backend_config.matrix_backend_override.is_some());
    }

    #[test]
    fn test_normalize_reactor_physics_task_map_coerces_solver_counters_for_rst_parser() {
        let mut document = DocumentMap::new();
        document.insert(
            "solver_settings".to_string(),
            HashMap::from([
                (
                    "max_iterations".to_string(),
                    Some(vec![Value::Integer(100)]),
                ),
                (
                    "refinement_steps".to_string(),
                    Some(vec![Value::Float(2.0)]),
                ),
            ]),
        );
        document.insert(
            "strategy_params".to_string(),
            HashMap::from([
                ("max_jac".to_string(), Some(vec![Value::Integer(3)])),
                ("max_damp_iter".to_string(), Some(vec![Value::Float(4.0)])),
            ]),
        );
        document.insert(
            "adaptive_strategy".to_string(),
            HashMap::from([
                ("version".to_string(), Some(vec![Value::Integer(1)])),
                (
                    "max_refinements".to_string(),
                    Some(vec![Value::Float(5.0)]),
                ),
            ]),
        );

        let normalized = normalize_reactor_physics_task_map(&document);

        assert!(matches!(
            normalized
                .get("solver_settings")
                .and_then(|section| section.get("max_iterations"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(100))
        ));
        assert!(matches!(
            normalized
                .get("solver_settings")
                .and_then(|section| section.get("refinement_steps"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(2))
        ));
        assert!(matches!(
            normalized
                .get("strategy_params")
                .and_then(|section| section.get("max_jac"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(3))
        ));
        assert!(matches!(
            normalized
                .get("strategy_params")
                .and_then(|section| section.get("max_damp_iter"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(4))
        ));
        assert!(matches!(
            normalized
                .get("adaptive_strategy")
                .and_then(|section| section.get("version"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(1))
        ));
        assert!(matches!(
            normalized
                .get("adaptive_strategy")
                .and_then(|section| section.get("max_refinements"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(5))
        ));
    }

    #[test]
    fn test_normalize_reactor_physics_task_map_unwraps_optional_string_counters() {
        let mut document = DocumentMap::new();
        document.insert(
            "solver_settings".to_string(),
            HashMap::from([(
                "max_iterations".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::String("120".to_string()))))]),
            )]),
        );
        document.insert(
            "adaptive_strategy".to_string(),
            HashMap::from([(
                "version".to_string(),
                Some(vec![Value::Optional(Some(Box::new(Value::String("2".to_string()))))]),
            )]),
        );

        let normalized = normalize_reactor_physics_task_map(&document);

        assert!(matches!(
            normalized
                .get("solver_settings")
                .and_then(|section| section.get("max_iterations"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(120))
        ));
        assert!(matches!(
            normalized
                .get("adaptive_strategy")
                .and_then(|section| section.get("version"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first()),
            Some(Value::Integer(2))
        ));
    }

    #[test]
    fn test_normalize_reactor_physics_task_map_handles_parsed_template_counters() {
        let mut parser = DocumentParser::new(task_content.to_string());
        parser
            .parse_document()
            .expect("task fixture should parse before normalization");
        let parsed = parser
            .get_result()
            .expect("parser should expose parsed task fixture");

        let normalized = normalize_reactor_physics_task_map(parsed);
        let value = normalized
            .get("solver_settings")
            .and_then(|section| section.get("max_iterations"))
            .and_then(|slot| slot.as_ref())
            .and_then(|values| values.first())
            .expect("normalized solver_settings.max_iterations should exist");

        match value {
            Value::Integer(value) => assert_eq!(*value, 100),
            other => panic!("max_iterations should normalize for the RST parser, got {other:?}"),
        }
    }

    #[test]
    fn test_normalized_template_solver_settings_parse_with_native_rst_parser() {
        let mut parser = DocumentParser::new(task_content.to_string());
        parser
            .parse_document()
            .expect("task fixture should parse before normalization");
        let parsed = parser
            .get_result()
            .expect("parser should expose parsed task fixture");

        let normalized = normalize_reactor_physics_task_map(parsed);
        let spec = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&normalized)
            .expect("normalized KiThe task fixture should satisfy the native RST parser");
        let options = spec
            .build_solver_options()
            .expect("normalized Lambdify settings should build native RST options");

        assert_eq!(spec.max_iterations, 100);
        assert!(matches!(
            options.generated_backend_config.backend_policy_override,
            Some(
                RustedSciThe::symbolic::codegen::codegen_backend_selection::BackendSelectionPolicy::LambdifyOnly
            )
        ));
        assert!(matches!(
            options.generated_backend_config.aot_build_policy,
            RustedSciThe::numerical::BVP_Damp::generated_solver_handoff::AotBuildPolicy::UseIfAvailable
        ));
        let solver_settings = normalized
            .get("solver_settings")
            .expect("normalized task should retain solver settings");
        for field_name in AOT_ONLY_SOLVER_FIELDS {
            assert!(
                !solver_settings.contains_key(*field_name),
                "Lambdify handoff must remove AOT-only field `{field_name}`"
            );
        }
    }

    #[test]
    fn test_lambdify_normalization_removes_aot_build_lifecycle() {
        let mut document = DocumentMap::new();
        document.insert(
            "solver_settings".to_string(),
            HashMap::from([
                (
                    "method".to_string(),
                    Some(vec![Value::String("Banded".to_string())]),
                ),
                (
                    "generated_backend".to_string(),
                    Some(vec![Value::String("banded_lambdify".to_string())]),
                ),
                (
                    "backend_policy".to_string(),
                    Some(vec![Value::String("lambdify_only".to_string())]),
                ),
                (
                    "aot_build_policy".to_string(),
                    Some(vec![Value::String("build_if_missing".to_string())]),
                ),
                (
                    "aot_c_compiler".to_string(),
                    Some(vec![Value::String("tcc".to_string())]),
                ),
            ]),
        );

        let normalized = normalize_reactor_physics_task_map(&document);
        let solver_settings = normalized
            .get("solver_settings")
            .expect("normalized task should retain solver settings");

        assert_eq!(
            solver_settings
                .get("backend_policy")
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first())
                .and_then(Value::as_string)
                .map(String::as_str),
            Some("lambdify_only")
        );
        assert!(!solver_settings.contains_key("aot_build_policy"));
        assert!(!solver_settings.contains_key("aot_c_compiler"));
    }

    #[test]
    fn test_normalize_reactor_physics_task_map_respects_adaptive_none_switch() {
        let mut document = DocumentMap::new();
        document.insert(
            "strategy_params".to_string(),
            HashMap::from([(
                "adaptive".to_string(),
                Some(vec![Value::Optional(None)]),
            )]),
        );
        document.insert(
            "adaptive_strategy".to_string(),
            HashMap::from([
                ("version".to_string(), Some(vec![Value::Integer(1)])),
                (
                    "max_refinements".to_string(),
                    Some(vec![Value::Integer(3)]),
                ),
            ]),
        );
        document.insert(
            "grid_refinement".to_string(),
            HashMap::from([(
                "pearson".to_string(),
                Some(vec![Value::Vector(vec![0.1, 1.5])]),
            )]),
        );

        let normalized = normalize_reactor_physics_task_map(&document);

        assert!(
            normalized.get("adaptive_strategy").is_none(),
            "strategy_params.adaptive: None must disable the native adaptive_strategy contract"
        );
        assert!(
            normalized.get("grid_refinement").is_none(),
            "grid_refinement is meaningful only when adaptive_strategy is enabled"
        );
    }

    #[test]
    fn test_adaptive_none_survives_gui_created_empty_sections_for_native_rst_parser() {
        let mut document = DocumentMap::new();
        document.insert(
            "solver_settings".to_string(),
            HashMap::from([
                (
                    "scheme".to_string(),
                    Some(vec![Value::String("forward".to_string())]),
                ),
                (
                    "method".to_string(),
                    Some(vec![Value::String("Banded".to_string())]),
                ),
                (
                    "strategy".to_string(),
                    Some(vec![Value::String("Damped".to_string())]),
                ),
                ("linear_sys_method".to_string(), Some(vec![Value::Optional(None)])),
                ("abs_tolerance".to_string(), Some(vec![Value::Float(1e-5)])),
                ("max_iterations".to_string(), Some(vec![Value::Usize(100)])),
                (
                    "loglevel".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::String(
                        "info".to_string(),
                    ))))]),
                ),
                ("dont_save_log".to_string(), Some(vec![Value::Boolean(true)])),
            ]),
        );
        document.insert(
            "bounds".to_string(),
            HashMap::from([
                (
                    "C".to_string(),
                    Some(vec![Value::Float(-10.0), Value::Float(10.0)]),
                ),
                (
                    "J".to_string(),
                    Some(vec![Value::Float(-1e20), Value::Float(1e20)]),
                ),
                (
                    "Teta".to_string(),
                    Some(vec![Value::Float(-100.0), Value::Float(100.0)]),
                ),
                (
                    "q".to_string(),
                    Some(vec![Value::Float(-1e20), Value::Float(1e20)]),
                ),
            ]),
        );
        document.insert(
            "rel_tolerance".to_string(),
            HashMap::from([
                ("C".to_string(), Some(vec![Value::Float(1e-5)])),
                ("J".to_string(), Some(vec![Value::Float(1e-5)])),
                ("Teta".to_string(), Some(vec![Value::Float(1e-5)])),
                ("q".to_string(), Some(vec![Value::Float(1e-5)])),
            ]),
        );
        document.insert(
            "strategy_params".to_string(),
            HashMap::from([
                (
                    "max_jac".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Integer(3))))]),
                ),
                (
                    "max_damp_iter".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Integer(10))))]),
                ),
                (
                    "damp_factor".to_string(),
                    Some(vec![Value::Optional(Some(Box::new(Value::Float(0.5))))]),
                ),
                ("adaptive".to_string(), Some(vec![Value::Optional(None)])),
            ]),
        );
        // Structured GUI rendering can create these sections even when the task
        // uses the legacy-but-valid `adaptive: None` switch.
        document.insert("adaptive_strategy".to_string(), HashMap::new());
        document.insert("grid_refinement".to_string(), HashMap::new());

        let normalized = normalize_reactor_physics_task_map(&document);
        let spec = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&normalized)
            .expect("adaptive: None should disable grid adaptation without requiring version");

        assert_eq!(spec.max_iterations, 100);
        assert_eq!(
            spec.strategy_params
                .expect("strategy params should be parsed")
                .adaptive,
            None
        );
    }

    #[test]
    fn test_parse_tolerance_and_bounds_alias_matches_legacy_spelling() {
        let mut reactor = SimpleReactorTask::new();

        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test_bounds_alias.txt");

        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        let mut legacy_parser = reactor
            .parse_file(Some(file_path.clone()))
            .expect("Failed to parse file");
        let _ = legacy_parser.parse_document();

        let mut alias_parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = alias_parser.parse_document();

        let legacy = reactor
            .parse_toleranse_and_bounds(&mut legacy_parser)
            .expect("Legacy spelling should still work");
        let alias = reactor
            .parse_tolerance_and_bounds(&mut alias_parser)
            .expect("New spelling should work too");

        assert_eq!(legacy.0.C, alias.0.C);
        assert_eq!(legacy.0.J, alias.0.J);
        assert_eq!(legacy.0.Teta, alias.0.Teta);
        assert_eq!(legacy.0.q, alias.0.q);
        assert_eq!(legacy.1.C, alias.1.C);
        assert_eq!(legacy.1.J, alias.1.J);
        assert_eq!(legacy.1.Teta, alias.1.Teta);
        assert_eq!(legacy.1.q, alias.1.q);
    }

    #[test]
    fn test_parse_grid_refinement_alias_matches_legacy_spelling() {
        let mut reactor = SimpleReactorTask::new();

        let legacy_task = task_content.to_string();
        let alias_task = task_content
            .replace("grcarsmooke", "grcar_smooke")
            .replace("doubleoints", "double_points");

        let mut legacy_parser = DocumentParser::new(legacy_task);
        legacy_parser
            .parse_document()
            .expect("Legacy spelling should parse");
        let mut alias_parser = DocumentParser::new(alias_task);
        alias_parser
            .parse_document()
            .expect("Corrected spelling should parse");

        let legacy = reactor.parse_parameters_with_exact_names(&mut legacy_parser);
        let alias = reactor.parse_parameters_with_exact_names(&mut alias_parser);

        assert!(legacy.is_ok());
        assert!(alias.is_ok());
    }

    #[test]
    fn test_parse_toleranse_and_bounds_accepts_scalar_pairs() {
        let mut document = DocumentMap::new();
        document.insert(
            "solver_settings".to_string(),
            HashMap::from([
                (
                    "scheme".to_string(),
                    Some(vec![Value::String("forward".to_string())]),
                ),
                (
                    "strategy".to_string(),
                    Some(vec![Value::String("Damped".to_string())]),
                ),
                (
                    "method".to_string(),
                    Some(vec![Value::String("Banded".to_string())]),
                ),
                ("abs_tolerance".to_string(), Some(vec![Value::Float(1e-7)])),
                ("max_iterations".to_string(), Some(vec![Value::Integer(100)])),
            ]),
        );
        document.insert(
            "bounds".to_string(),
            HashMap::from([
                (
                    "C".to_string(),
                    Some(vec![Value::Float(-10.0), Value::Float(10.0)]),
                ),
                (
                    "J".to_string(),
                    Some(vec![Value::Float(-1e20), Value::Float(1e20)]),
                ),
                (
                    "Teta".to_string(),
                    Some(vec![Value::Float(-100.0), Value::Float(100.0)]),
                ),
                (
                    "q".to_string(),
                    Some(vec![Value::Float(-1e20), Value::Float(1e20)]),
                ),
            ]),
        );
        document.insert(
            "rel_tolerance".to_string(),
            HashMap::from([
                ("C".to_string(), Some(vec![Value::Float(1e-7)])),
                ("J".to_string(), Some(vec![Value::Float(1e-7)])),
                ("Teta".to_string(), Some(vec![Value::Float(1e-7)])),
                ("q".to_string(), Some(vec![Value::Float(1e-7)])),
            ]),
        );

        let spec = task_parser_damped::parse_bvp_damped_solver_settings_from_document(&document)
            .expect("Failed to parse native damped solver settings");

        assert_eq!(spec.abs_tolerance, 1e-7);
        assert_eq!(
            spec.rel_tolerance.as_ref().and_then(|map| map.get("C")).copied(),
            Some(1e-7)
        );
        assert_eq!(
            spec.bounds.as_ref().and_then(|map| map.get("C")).copied(),
            Some((-10.0, 10.0))
        );
    }

    #[test]
    fn test_validate_solver_output_matrix_rejects_non_finite_values() {
        let matrix = DMatrix::from_vec(2, 2, vec![1.0, 2.0, 3.0, f64::NAN]);
        let result = validate_solver_output_matrix(&matrix);
        assert!(result.is_err());
        match result {
            Err(ReactorError::InvalidNumericValue(msg)) => {
                assert!(msg.contains("non-finite value"));
            }
            _ => panic!("Expected InvalidNumericValue error"),
        }
    }

    #[test]
    fn test_parse_parameters_rejects_non_finite_pressure() {
        let mut reactor = SimpleReactorTask::new();
        let task = task_content.replace("P: 1e6", "P: NaN");
        let mut parser = DocumentParser::new(task);
        let result = reactor.parse_parameters_with_exact_names(&mut parser);
        assert!(matches!(result, Err(ReactorError::InvalidNumericValue(_))));
    }

    #[test]
    fn test_parse_parameters_with_exact_names_is_idempotent_for_preparsed_document() {
        let mut reactor = SimpleReactorTask::new();
        let mut parser = DocumentParser::new(task_content.to_string());
        parser.parse_document().expect("Failed to parse document");

        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("First parse should succeed");
        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("Second parse should also succeed");

        assert_eq!(reactor.kindata.substances, vec!["HMX", "HMXprod"]);
        assert_eq!(reactor.kindata.vec_of_equations, vec!["HMX=>10HMXprod"]);
        assert_eq!(reactor.boundary_condition.get("HMX"), Some(&0.999));
        assert_eq!(reactor.boundary_condition.get("HMXprod"), Some(&0.001));
    }

    #[test]
    fn test_setup_bvp_from_file() {
        let mut reactor = SimpleReactorTask::new();
        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");
        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        // println!("parser {:?}", parser);

        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("Failed to parse parameters");

        reactor.setup_bvp().expect("Failed to setup BVP");
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            info!("reaction {} rate {}", eq, rate);
        }
        info!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            info!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        info!("bc {:?}", bc);
        info!(" unknowns{:?}", reactor.solver.unknowns);
    }

    #[test]
    fn test_solve_from_map() {
        let mut reactor = SimpleReactorTask::new();
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");
        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        reactor
            .solve_from_map(parser)
            .expect("Failed to solve from map");
        // println!("parser {:?}", parser);
    }
}
