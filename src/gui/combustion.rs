//! # Combustion GUI Module
//!
//! This module provides a comprehensive graphical user interface for configuring and running
//! chemical reactor boundary value problems (BVPs), particularly gas-phase combustion and
//! steady-state plug-flow reactor simulations.
//!
//! ## Architecture Overview
//!
//! The GUI is built using the egui framework and provides:
//! - **Interactive document editing**: Users can modify reactor parameters through a structured interface
//! - **File I/O operations**: Load, save, and save-as functionality for task configuration files
//! - **Problem type selection**: Support for different reactor problem variants
//! - **Real-time validation**: Input validation and type conversion for numerical parameters
//! - **Integrated help system**: Context-sensitive help based on selected problem type
//! - **Calculation execution**: Direct integration with BVP solvers
//!
//! ## Data Flow
//!
//! 1. **Configuration Loading**: Task files are parsed into DocumentMap structure
//! 2. **Interactive Editing**: Users modify parameters through GUI widgets
//! 3. **Validation & Conversion**: Input values are validated and converted to appropriate types
//! 4. **Calculation Execution**: Modified parameters are passed to appropriate solvers
//! 5. **Results Processing**: Solver results are processed and displayed
//!
//! ## Key Components
//!
//! - `CombustionApp`: Main application state and UI logic
//! - `ProblemsEnum`: Problem type enumeration for different reactor configurations
//! - `render_section()`: Dynamic UI generation for document sections
//! - `render_value()`: Type-specific widget rendering for different value types
//!
//! ## Supported Problem Types
//!
//! - **BVPSimple**: Gas-phase combustion and plug-flow reactor BVP problems
//! - **None**: Empty configuration for custom problem setup

use crate::ReactorsBVP::SimpleReactorBVP::SimpleReactorTask;
use crate::ReactorsBVP::task_parser_reactor_BVP::{
    SIMPLE_BVP_TEMPLATE, validate_reactor_bvp_task_map,
};
use crate::ReactorsBVP::task_value_conversion::{UsizeConversionError, try_value_as_usize};
use crate::cli::reactor_help::REACTOR_ENG_HELPER;
use crate::gui::bvp_gui_config::{
    AotToolchainRequest, BvpGuiConfig, BvpGuiMigrationReport, GuiAotCodegenBackend,
    GuiExecutionBackend, aot_toolchain_request, normalize_bvp_solver_backend_section,
    validate_aot_toolchain_fields,
};
pub use crate::gui::document_lifecycle::{
    DocumentLifecycleState, GuiFileOperationKind, GuiFileOperationResult,
};
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser, Value};
use eframe::egui;
use log::{info, warn};
use std::collections::HashMap;
use std::sync::atomic::AtomicBool;
use std::sync::{Arc, mpsc};
use tabled::{Table, Tabled};

#[path = "combustion_execution.rs"]
mod combustion_execution;
#[path = "combustion_lifecycle.rs"]
mod combustion_lifecycle;
/*
this is value enum from RustedSciThe::command_interpreter::task_parser
pub enum Value {
    String(String),
    Float(f64),
    Integer(i64),
    Usize(usize),
    Vector(Vec<f64>),
    Boolean(bool),
    Optional(Option<Box<Value>>),
}
pub type DocumentMap = HashMap<String, SectionMap>;
type SectionMap = HashMap<String, Option<Vec<Value>>>;
DocumentParser is a struct that contains result of parsing and
 has a method get_result() -> Option<DocumentMap>
*/

/// Enumeration of supported reactor problem types.
///
/// This enum defines the different types of reactor problems that can be configured
/// and solved through the GUI. Each variant corresponds to a specific mathematical
/// model and solver configuration.
///
/// # Variants
///
/// * `None` - Empty configuration, allows users to build custom problems from scratch
/// * `BVPSimple` - Gas-phase combustion and steady-state plug-flow reactor problems
///   using boundary value problem formulation with dimensionless variables
#[derive(Default, Debug, Clone, PartialEq)]
pub enum ProblemsEnum {
    #[default]
    None,
    BVPSimple,
}
impl ProblemsEnum {
    /// Returns the template configuration string for the selected problem type.
    ///
    /// Each problem type has an associated template that provides default values
    /// and structure for the configuration parameters. This template is parsed
    /// into a DocumentMap for GUI editing.
    ///
    /// # Returns
    ///
    /// * `Some(&str)` - Template string for the problem type
    /// * `None` - No template available (empty configuration)
    pub fn get_problem(&self) -> Option<&str> {
        match self {
            ProblemsEnum::None => None,
            ProblemsEnum::BVPSimple => Some(SIMPLE_BVP_TEMPLATE),
        }
    }
}
impl std::fmt::Display for ProblemsEnum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            ProblemsEnum::None => "None",
            ProblemsEnum::BVPSimple => "BVP Simple",
        };
        write!(f, "{}", s)
    }
}

/// Renders an interactive GUI widget for a specific Value type.
///
/// This function provides type-specific rendering for different Value variants,
/// creating appropriate egui widgets for user interaction. Each widget type
/// is optimized for the specific data type it represents.
///
/// # Arguments
///
/// * `ui` - The egui UI context for rendering widgets
/// * `value` - Mutable reference to the Value being edited
/// * `field_name` - Name of the field for logging and identification
///
/// # Widget Types
///
/// * `String` - Single-line text input
/// * `Float` - Drag value widget with 0.1 step size
/// * `Integer` - Drag value widget with 1 step size
/// * `Usize` - Non-negative drag value widget
/// * `Vector` - Multiple drag widgets for each vector element
/// * `Boolean` - Checkbox widget
/// * `Optional` - Read/edit an existing typed payload without inventing its type
fn render_value(ui: &mut egui::Ui, value: &mut Value, field_name: &str) {
    match value {
        Value::String(s) => {
            if ui.text_edit_singleline(s).changed() {
                info!("Field '{}' changed to: {}", field_name, s);
            }
        }
        Value::Float(f) => {
            if ui
                .add(egui::DragValue::new(f).speed(0.1).custom_formatter(|n, _| {
                    if n.abs() < 1e-3 && n != 0.0 {
                        format!("{:.2e}", n)
                    } else if n.abs() >= 100.0 {
                        format!("{:.1}", n)
                    } else {
                        format!("{:.6}", n)
                    }
                }))
                .changed()
            {
                info!("Field '{}' changed to: {}", field_name, f);
            }
        }
        Value::Integer(i) => {
            if ui.add(egui::DragValue::new(i).speed(1)).changed() {
                info!("Field '{}' changed to: {}", field_name, i);
            }
        }
        Value::Usize(u) => {
            if ui.add(egui::DragValue::new(u).speed(1)).changed() {
                info!("Field '{}' changed to: {}", field_name, u);
            }
        }
        Value::Vector(vec) => {
            for (i, v) in vec.iter_mut().enumerate() {
                if ui
                    .add(egui::DragValue::new(v).speed(0.1).custom_formatter(|n, _| {
                        if n.abs() < 1e-3 && n != 0.0 {
                            format!("{:.2e}", n)
                        } else if n.abs() >= 100.0 {
                            format!("{:.1}", n)
                        } else {
                            format!("{:.6}", n)
                        }
                    }))
                    .changed()
                {
                    info!("Field '{}[{}]' changed to: {}", field_name, i, v);
                }
            }
        }
        Value::Boolean(b) => {
            if ui.checkbox(b, "").changed() {
                info!("Field '{}' changed to: {}", field_name, b);
            }
        }
        Value::Optional(opt) => {
            ui.horizontal(|ui| {
                if let Some(inner) = opt {
                    ui.label("Some:");
                    render_value(ui, inner, field_name);
                } else {
                    // A generic `None` carries no information about its inner
                    // type. Typed optional controls are responsible for
                    // constructing parser-valid `Some` values.
                    ui.label("None");
                }
            });
        }
    }
}

/// Renders a collapsible section containing multiple fields with add/delete functionality.
///
/// This function creates a collapsible header containing all fields in a section,
/// with individual delete buttons for each field and controls for adding new fields.
/// It handles the complex state management required for safe field deletion during iteration.
///
/// # Arguments
///
/// * `ui` - The egui UI context for rendering
/// * `name` - Display name for the section header
/// * `section` - Mutable reference to the section's field data
/// * `new_field_name` - Buffer for new field name input
/// * `new_field_value` - Buffer for new field value input
///
/// # Returns
///
/// A tuple containing:
/// * `bool` - Whether a new field should be added
/// * `bool` - Whether the entire section should be deleted
/// * `Vec<String>` - Field names whose deletion was requested
///
/// # Safety
///
/// Field deletion is handled safely by collecting fields to delete first,
/// then removing them outside the iteration to avoid borrowing conflicts.
fn render_section(
    ui: &mut egui::Ui,
    name: &str,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    new_field_name: &mut String,
    new_field_value: &mut String,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) -> (bool, bool, Vec<String>) {
    let mut add_field = false;
    let mut delete_section = false;
    let mut requested_field_deletions = Vec::new();

    ui.push_id(name, |ui| {
        ui.horizontal(|ui| {
            let section_header =
                egui::CollapsingHeader::new(name)
                    .default_open(false)
                    .show(ui, |ui| {
                        let mut fields_to_delete = Vec::new();
                        let mut field_names = sorted_section_field_names(section);
                        for field_name in field_names.drain(..) {
                            let Some(maybe_values) = section.get_mut(&field_name) else {
                                continue;
                            };
                            ui.push_id((name, field_name.as_str()), |ui| {
                                ui.horizontal(|ui| {
                                    let field_label = ui.label(&field_name);

                                    // Update help info on hover.
                                    if field_label.hovered() {
                                        if let Some(map) = help_map {
                                            let help_key = format!("{}.{}", name, field_name);
                                            if let Some((_, help_text)) =
                                                map.iter().find(|(key, _)| *key == help_key)
                                            {
                                                *current_help = help_text.to_string();
                                            }
                                        }
                                    }

                                    if let Some(values) = maybe_values {
                                        // Special handling for Vec<String> fields like substances.
                                        if values.len() > 1
                                            && values.iter().all(|v| matches!(v, Value::String(_)))
                                        {
                                            let mut combined_text = values
                                                .iter()
                                                .filter_map(|v| {
                                                    if let Value::String(s) = v {
                                                        Some(s.as_str())
                                                    } else {
                                                        None
                                                    }
                                                })
                                                .collect::<Vec<_>>()
                                                .join(", ");

                                            if ui.text_edit_singleline(&mut combined_text).changed()
                                            {
                                                let new_substances: Vec<Value> = combined_text
                                                    .split(',')
                                                    .map(|s| Value::String(s.trim().to_string()))
                                                    .filter(|v| {
                                                        if let Value::String(s) = v {
                                                            !s.is_empty()
                                                        } else {
                                                            true
                                                        }
                                                    })
                                                    .collect();
                                                *values = new_substances;
                                                info!(
                                                    "Field '{}' changed to: {}",
                                                    field_name, combined_text
                                                );
                                            }
                                        } else {
                                            for v in values.iter_mut() {
                                                render_value(ui, v, &field_name);
                                            }
                                        }
                                    } else {
                                        ui.label("None");
                                    }

                                    if ui.button("❌").clicked() {
                                        fields_to_delete.push(field_name.clone());
                                    }
                                });
                            });
                        }

                        requested_field_deletions.extend(fields_to_delete);

                        ui.separator();
                        ui.horizontal(|ui| {
                            ui.label("New field:");
                            ui.text_edit_singleline(new_field_name);
                            ui.label("Value:");
                            ui.text_edit_singleline(new_field_value);
                            if ui.button("Add Field").clicked() {
                                add_field = true;
                            }
                        });
                    });

            // Update help info on section header hover
            if section_header.header_response.hovered() {
                if let Some(map) = help_map {
                    if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == name) {
                        *current_help = help_text.to_string();
                    }
                }
            }

            if ui.button("❌").clicked() {
                delete_section = true;
            }
        });
    });

    (add_field, delete_section, requested_field_deletions)
}

/// Returns section field names in a deterministic order for raw document rendering.
pub(crate) fn sorted_section_field_names(
    section: &HashMap<String, Option<Vec<Value>>>,
) -> Vec<String> {
    let mut field_names = section.keys().cloned().collect::<Vec<_>>();
    field_names.sort();
    field_names
}

/// Returns document section names in a deterministic order for raw GUI rendering.
pub(crate) fn sorted_document_section_names(document: &DocumentMap) -> Vec<String> {
    let mut section_names = document.keys().cloned().collect::<Vec<_>>();
    section_names.sort();
    section_names
}

const BVP_GRID_REFINEMENT_METHODS: &[&str] =
    &["pearson", "grcarsmooke", "twopnt", "easy", "doubleoints"];

/// Canonicalizes grid-refinement method names accepted by older task files.
///
/// RustedSciThe currently parses `grcarsmooke` and `doubleoints`; KiThe also
/// accepts the more readable aliases in the GUI and rewrites them at the
/// structured-editor boundary.
fn canonical_grid_refinement_method_name(raw: &str) -> Option<&'static str> {
    let normalized = raw.trim().to_ascii_lowercase().replace(['-', ' '], "_");
    match normalized.as_str() {
        "pearson" => Some("pearson"),
        "grcarsmooke" | "grcar_smooke" => Some("grcarsmooke"),
        "twopnt" | "two_point" | "two_points" => Some("twopnt"),
        "easy" => Some("easy"),
        "doubleoints" | "double_points" | "doublepoints" => Some("doubleoints"),
        _ => None,
    }
}

fn grid_refinement_method_label(method: &str) -> &'static str {
    match method {
        "pearson" => "Pearson",
        "grcarsmooke" => "Grcar-Smooke",
        "twopnt" => "Two-point",
        "easy" => "Easy",
        "doubleoints" => "Double points",
        _ => "Pearson",
    }
}

fn default_grid_refinement_params(method: &str) -> Vec<f64> {
    match method {
        "pearson" => vec![0.1, 1.5],
        "grcarsmooke" => vec![0.05, 0.05, 1.25],
        "twopnt" => vec![0.05, 0.05, 1.25],
        "easy" => vec![0.05],
        "doubleoints" => Vec::new(),
        _ => vec![0.1, 1.5],
    }
}

fn grid_refinement_slot_to_params(slot: Option<&Vec<Value>>) -> Option<Vec<f64>> {
    let values = slot?;
    if let Some(vector) = values.first().and_then(|value| value.as_vector()) {
        return Some(vector.clone());
    }

    let mut params = Vec::with_capacity(values.len());
    for value in values {
        params.push(value.as_float()?);
    }
    Some(params)
}

/// Reshapes a grid-refinement section into the single-method contract used by RST.
fn normalize_grid_refinement_section(section: &mut HashMap<String, Option<Vec<Value>>>) -> String {
    let mut selected_method = None;
    let mut selected_params = None;

    for method in BVP_GRID_REFINEMENT_METHODS {
        if let Some((key, slot)) = section
            .iter()
            .find(|(key, _)| canonical_grid_refinement_method_name(key) == Some(*method))
        {
            selected_method = Some(*method);
            selected_params = grid_refinement_slot_to_params(slot.as_ref());
            if key.as_str() != *method {
                info!("Grid refinement alias '{}' normalized to '{}'", key, method);
            }
            break;
        }
    }

    let method = selected_method.unwrap_or("pearson");
    let mut params = default_grid_refinement_params(method);
    if let Some(existing_params) = selected_params {
        for (target, source) in params.iter_mut().zip(existing_params.into_iter()) {
            if source.is_finite() {
                *target = source;
            }
        }
    }

    section.clear();
    section.insert(method.to_string(), Some(vec![Value::Vector(params)]));
    method.to_string()
}

/// Returns whether a compatibility value explicitly disables adaptive refinement.
#[allow(dead_code)]
fn value_is_explicit_none(value: &Value) -> bool {
    matches!(value, Value::Optional(None))
        || matches!(value, Value::String(value) if value.trim().eq_ignore_ascii_case("none"))
}

/// Reads and removes the obsolete `strategy_params.adaptive` GUI switch.
///
/// RustedSciThe enables adaptive refinement from the presence of the
/// `adaptive_strategy` section, not from this field. Older task documents may
/// still carry `adaptive: None`, so the marker is consumed once and translated
/// into the native section-presence contract.
#[allow(dead_code)]
fn take_legacy_adaptive_marker(document: &mut DocumentMap) -> Option<bool> {
    let strategy_params = document.get_mut("strategy_params")?;
    let values = strategy_params.remove("adaptive")??;
    let value = values.first()?;
    Some(!value_is_explicit_none(value))
}

/// Enables or disables the complete native adaptive-grid contract.
///
/// Version 1 is currently the only production strategy exposed by
/// RustedSciThe. Enabling adaptation therefore seeds it automatically together
/// with a bounded refinement count and one valid grid-refinement method.
fn set_bvp_adaptive_grid_enabled(document: &mut DocumentMap, enabled: bool) {
    if let Some(strategy_params) = document.get_mut("strategy_params") {
        strategy_params.remove("adaptive");
    }

    if !enabled {
        document.remove("adaptive_strategy");
        document.remove("grid_refinement");
        return;
    }

    let adaptive_strategy = ensure_section_mut(document, "adaptive_strategy");
    adaptive_strategy.insert("version".to_string(), Some(vec![Value::Usize(1)]));
    let _ = ensure_field_slot(adaptive_strategy, "max_refinements", Value::Usize(3));
    let _ = normalize_usize_field(adaptive_strategy, "max_refinements");

    let grid_refinement = ensure_section_mut(document, "grid_refinement");
    normalize_grid_refinement_section(grid_refinement);
}

/// Canonicalizes adaptive-grid state loaded from old or partially edited tasks.
#[allow(dead_code)]
fn normalize_bvp_adaptive_grid_settings(document: &mut DocumentMap) {
    let compatibility_marker = take_legacy_adaptive_marker(document);
    let has_native_payload = document
        .get("adaptive_strategy")
        .is_some_and(|section| !section.is_empty())
        || document
            .get("grid_refinement")
            .is_some_and(|section| !section.is_empty());
    let enabled = compatibility_marker.unwrap_or(has_native_payload);
    set_bvp_adaptive_grid_enabled(document, enabled);
}

/// Renders the one coherent adaptive-grid control surface accepted by RST.
fn render_bvp_adaptive_grid_panel(
    ui: &mut egui::Ui,
    document: &mut DocumentMap,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    let mut enabled = document.contains_key("adaptive_strategy");

    ui.group(|ui| {
        ui.vertical(|ui| {
            let response = ui.checkbox(&mut enabled, "Adaptive grid refinement");
            if response.hovered() {
                if let Some(map) = help_map {
                    if let Some((_, help_text)) =
                        map.iter().find(|(key, _)| *key == "adaptive_strategy")
                    {
                        *current_help = help_text.to_string();
                    }
                }
            }
            if response.changed() {
                set_bvp_adaptive_grid_enabled(document, enabled);
                info!("Adaptive grid refinement enabled: {}", enabled);
            }

            if enabled {
                ui.label("Adaptive strategy version: 1");
                egui::Grid::new("bvp_adaptive_strategy_grid")
                    .num_columns(2)
                    .striped(true)
                    .show(ui, |ui| {
                        let adaptive_strategy = ensure_section_mut(document, "adaptive_strategy");
                        let current_result = adaptive_strategy
                            .get("max_refinements")
                            .and_then(|slot| slot.as_ref())
                            .and_then(|values| values.first())
                            .map(try_value_as_usize)
                            .unwrap_or(Ok(3));
                        ui.label("Maximum grid refinements");
                        match current_result {
                            Ok(mut max_refinements) => {
                                if ui
                                    .add(egui::DragValue::new(&mut max_refinements).speed(1))
                                    .changed()
                                {
                                    adaptive_strategy.insert(
                                        "max_refinements".to_string(),
                                        Some(vec![Value::Usize(max_refinements)]),
                                    );
                                    info!(
                                        "Maximum grid refinements changed to: {}",
                                        max_refinements
                                    );
                                }
                            }
                            Err(error) => {
                                ui.horizontal(|ui| {
                                    ui.colored_label(egui::Color32::RED, error.to_string());
                                    if ui.button("Reset to 3").clicked() {
                                        adaptive_strategy.insert(
                                            "max_refinements".to_string(),
                                            Some(vec![Value::Usize(3)]),
                                        );
                                    }
                                });
                            }
                        }
                        ui.end_row();
                    });

                let (_, delete_requested) = render_grid_refinement_section(
                    ui,
                    ensure_section_mut(document, "grid_refinement"),
                    help_map,
                    current_help,
                );
                if delete_requested {
                    enabled = false;
                    set_bvp_adaptive_grid_enabled(document, false);
                }
            }
        });
    });
}

/// Renders grid refinement as a finite method selector plus typed parameters.
///
/// A raw map is too easy to make invalid here because RustedSciThe accepts only
/// one active method. This editor keeps the document in the exact shape the
/// backend parser expects.
fn render_grid_refinement_section(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) -> (bool, bool) {
    let mut delete_section = false;
    let current_method = normalize_grid_refinement_section(section);

    ui.horizontal(|ui| {
        let section_header = egui::CollapsingHeader::new("grid_refinement")
            .default_open(true)
            .show(ui, |ui| {
                egui::Grid::new("bvp_grid_refinement_strategy_grid")
                    .num_columns(2)
                    .striped(true)
                    .show(ui, |ui| {
                        let label_response = ui.label("Grid refinement strategy");
                        if label_response.hovered() {
                            if let Some(map) = help_map {
                                if let Some((_, help_text)) =
                                    map.iter().find(|(key, _)| *key == "grid_refinement")
                                {
                                    *current_help = help_text.to_string();
                                }
                            }
                        }

                        let mut selected_method = current_method.clone();
                        egui::ComboBox::from_id_salt("grid_refinement.method")
                            .selected_text(grid_refinement_method_label(&selected_method))
                            .show_ui(ui, |ui| {
                                for method in BVP_GRID_REFINEMENT_METHODS {
                                    ui.selectable_value(
                                        &mut selected_method,
                                        (*method).to_string(),
                                        grid_refinement_method_label(method),
                                    );
                                }
                            });
                        ui.end_row();

                        if selected_method != current_method {
                            section.clear();
                            section.insert(
                                selected_method.clone(),
                                Some(vec![Value::Vector(default_grid_refinement_params(
                                    &selected_method,
                                ))]),
                            );
                            info!("Grid refinement strategy changed to: {}", selected_method);
                        }

                        let params = section
                            .get_mut(&selected_method)
                            .and_then(|slot| slot.as_mut())
                            .and_then(|values| values.first_mut())
                            .and_then(|value| match value {
                                Value::Vector(params) => Some(params),
                                _ => None,
                            })
                            .expect("grid refinement section should contain a vector payload");

                        ui.label("Parameters");
                        if params.is_empty() {
                            ui.label("No parameters");
                        } else {
                            ui.horizontal(|ui| {
                                for (index, value) in params.iter_mut().enumerate() {
                                    if ui
                                        .add(
                                            egui::DragValue::new(value)
                                                .speed(0.01)
                                                .prefix(format!("p{}=", index + 1)),
                                        )
                                        .changed()
                                    {
                                        info!(
                                            "Grid refinement parameter {} changed to: {}",
                                            index + 1,
                                            value
                                        );
                                    }
                                }
                            });
                        }
                        ui.end_row();
                    });
            });

        if section_header.header_response.hovered() {
            if let Some(map) = help_map {
                if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == "grid_refinement")
                {
                    *current_help = help_text.to_string();
                }
            }
        }

        if ui.button("Delete").clicked() {
            delete_section = true;
        }
    });

    (false, delete_section)
}

/// Reads the current grid-refinement selection without mutating the document.
#[allow(dead_code)]
fn current_grid_refinement_state(
    section: &HashMap<String, Option<Vec<Value>>>,
) -> (String, Vec<f64>) {
    for method in BVP_GRID_REFINEMENT_METHODS {
        if let Some((_key, slot)) = section
            .iter()
            .find(|(key, _)| canonical_grid_refinement_method_name(key) == Some(*method))
        {
            let params = grid_refinement_slot_to_params(slot.as_ref())
                .unwrap_or_else(|| default_grid_refinement_params(method));
            return ((*method).to_string(), params);
        }
    }

    let method = "pearson".to_string();
    (method.clone(), default_grid_refinement_params(&method))
}

/// Serializes a parser value using the same DSL that `DocumentParser` accepts.
fn value_to_document_string(value: &Value) -> String {
    match value {
        Value::String(value) => value.clone(),
        Value::Float(value) => value.to_string(),
        Value::Integer(value) => value.to_string(),
        Value::Usize(value) => value.to_string(),
        Value::Boolean(value) => value.to_string(),
        Value::Vector(values) => {
            let values = values
                .iter()
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            format!("[{}]", values)
        }
        Value::Optional(Some(inner)) => format!("Some({})", value_to_document_string(inner)),
        Value::Optional(None) => "None".to_string(),
    }
}

/// Serializes all values assigned to one document key as one parser-friendly line.
///
/// `DocumentParser` accepts `key: a, b`, but rejects repeating the same key
/// several times in one section. The GUI model may still store scalar pairs as
/// `Some(vec![Value::Float(a), Value::Float(b)])`, so saving must collapse them
/// back into a single comma-separated DSL value.
fn values_to_document_string(values: &[Value]) -> String {
    values
        .iter()
        .map(value_to_document_string)
        .collect::<Vec<_>>()
        .join(", ")
}

/// Ensures that a section exists and returns a mutable reference to it.
///
/// The BVP solver panel seeds missing solver fields on demand, but it never
/// overwrites existing user edits.
fn ensure_section_mut<'a>(
    document: &'a mut DocumentMap,
    section_name: &str,
) -> &'a mut HashMap<String, Option<Vec<Value>>> {
    document
        .entry(section_name.to_string())
        .or_insert_with(HashMap::new)
}

/// Insert a user-created section without replacing existing task data.
pub(crate) fn insert_section_if_absent(
    document: &mut DocumentMap,
    section_name: String,
) -> Result<(), String> {
    match document.entry(section_name.clone()) {
        std::collections::hash_map::Entry::Vacant(entry) => {
            entry.insert(HashMap::new());
            Ok(())
        }
        std::collections::hash_map::Entry::Occupied(_) => Err(format!(
            "Section `{}` already exists; existing data was preserved",
            section_name
        )),
    }
}

/// Insert a user-created field without silently replacing its current value.
pub(crate) fn insert_field_if_absent(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: String,
    value: Value,
) -> Result<(), String> {
    match section.entry(field_name.clone()) {
        std::collections::hash_map::Entry::Vacant(entry) => {
            entry.insert(Some(vec![value]));
            Ok(())
        }
        std::collections::hash_map::Entry::Occupied(_) => Err(format!(
            "Field `{}` already exists; existing data was preserved",
            field_name
        )),
    }
}

/// One destructive document edit awaiting explicit user confirmation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum DocumentDeleteRequest {
    Section { section: String },
    Field { section: String, field: String },
}

impl DocumentDeleteRequest {
    fn description(&self) -> String {
        match self {
            Self::Section { section } => format!("section `{}`", section),
            Self::Field { section, field } => {
                format!("field `{}.{}`", section, field)
            }
        }
    }
}

/// Apply an already-authorized deletion and report whether data was removed.
pub(crate) fn apply_document_deletion(
    document: &mut DocumentMap,
    request: &DocumentDeleteRequest,
) -> bool {
    match request {
        DocumentDeleteRequest::Section { section } => document.remove(section).is_some(),
        DocumentDeleteRequest::Field { section, field } => document
            .get_mut(section)
            .is_some_and(|section_map| section_map.remove(field).is_some()),
    }
}

/// Ensures that a field slot exists and seeds a single default value if needed.
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

/// Seeds the canonical BVP solver backend defaults.
///
/// This keeps the structured GUI aligned with the lower RustedSciThe solver
/// API and gives users a visible starting point for both lambdify and AOT
/// routes.
fn ensure_bvp_solver_backend_defaults(document: &mut DocumentMap) {
    let section = ensure_section_mut(document, "solver_settings");
    let _ = normalize_matrix_backend_aliases(section);
    let defaults = [
        ("scheme", Value::String("forward".to_string())),
        ("strategy", Value::String("Damped".to_string())),
        ("linear_sys_method", Value::Optional(None)),
        (
            "generated_backend",
            Value::String("banded_lambdify".to_string()),
        ),
        ("symbolic_backend", Value::String("AtomView".to_string())),
        ("banded_linear_solver", Value::String("auto".to_string())),
        ("refinement_steps", Value::Usize(5)),
        ("abs_tolerance", Value::Float(1e-7)),
        ("max_iterations", Value::Usize(100)),
        (
            "loglevel",
            Value::Optional(Some(Box::new(Value::String("info".to_string())))),
        ),
        ("dont_save_log", Value::Boolean(true)),
    ];

    for (field_name, default_value) in defaults {
        let _ = ensure_field_slot(section, field_name, default_value);
    }

    let _ = normalize_usize_field(section, "max_iterations");
    let _ = normalize_usize_field(section, "refinement_steps");
    normalize_bvp_solver_backend_section(section);
}

/// Seeds the BVP document once during load or problem construction.
///
/// Rendering should stay read-only. This helper canonicalizes solver defaults
/// and the adaptive-grid sections before widgets inspect the document.
fn seed_bvp_document_defaults(document: &mut DocumentMap) {
    // Reuse the typed GUI snapshot as the single canonical source for the
    // seeded structured model, then write it back in one pass.
    let (config, _) = BvpGuiConfig::from_document(document);
    config.apply_to_document(document);
}

/// Normalizes numeric fields that must stay `usize` for the solver backend.
///
/// The GUI can load or temporarily hold numeric values in different scalar
/// variants, but the solver contract expects these fields to remain non-negative
/// integers.
fn normalize_usize_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) -> Result<(), UsizeConversionError> {
    let Some(values) = section.get_mut(field_name).and_then(|slot| slot.as_mut()) else {
        return Ok(());
    };

    if let Some(first) = values.first_mut() {
        let value = try_value_as_usize(first)?;
        *first = Value::Usize(value);
        values.truncate(1);
    }
    Ok(())
}

/// Canonicalizes every integer field crossing from the GUI into reactor code.
///
/// Missing fields remain the responsibility of the typed physics/solver
/// parsers. Present fields must satisfy the same non-lossy conversion contract
/// regardless of which editor or compatibility document produced them.
pub(crate) fn normalize_bvp_usize_fields(document: &mut DocumentMap) -> Result<(), String> {
    const FIELDS: &[(&str, &str)] = &[
        ("process_conditions", "n_steps"),
        ("solver_settings", "max_iterations"),
        ("solver_settings", "refinement_steps"),
        ("strategy_params", "max_jac"),
        ("strategy_params", "max_damp_iter"),
        ("adaptive_strategy", "version"),
        ("adaptive_strategy", "max_refinements"),
    ];

    for (section_name, field_name) in FIELDS {
        let Some(section) = document.get_mut(*section_name) else {
            continue;
        };
        normalize_usize_field(section, field_name).map_err(|error| {
            format!("invalid integer field `{section_name}.{field_name}`: {error}")
        })?;
    }
    Ok(())
}

fn string_field_value(
    section: &HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) -> Option<String> {
    section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .and_then(|value| value.as_string().cloned())
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

fn canonical_matrix_backend(value: &str) -> Option<&'static str> {
    match value.trim().to_ascii_lowercase().as_str() {
        "banded" => Some("Banded"),
        "sparse" => Some("Sparse"),
        _ => None,
    }
}

/// Keeps the two RST matrix-backend aliases synchronized behind one GUI field.
///
/// `method` is the required solver setting while `matrix_backend` is the
/// generated-backend override. They select the same Sparse/Banded storage route
/// for reactor BVP tasks, so exposing both independently only permits invalid
/// disagreement. The required `method` field wins when loading old documents.
fn normalize_matrix_backend_aliases(
    section: &mut HashMap<String, Option<Vec<Value>>>,
) -> &'static str {
    let method = string_field_value(section, "method")
        .as_deref()
        .and_then(canonical_matrix_backend);
    let matrix_backend = string_field_value(section, "matrix_backend")
        .as_deref()
        .and_then(canonical_matrix_backend);
    let backend = method.or(matrix_backend).unwrap_or("Banded");

    set_string_field(section, "method", backend);
    set_string_field(section, "matrix_backend", backend);
    backend
}

/// Executables that KiThe permits a task document to select directly.
///
fn optional_string_field_value(
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

fn normalize_optional_string_field(
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
) {
    let current_value = optional_string_field_value(section, field_name);
    let is_valid_shape = section
        .get(field_name)
        .and_then(|slot| slot.as_ref())
        .and_then(|values| values.first())
        .map(|value| match value {
            Value::Optional(None) => true,
            Value::Optional(Some(inner)) => inner.as_string().is_some(),
            Value::String(_) => true,
            _ => false,
        })
        .unwrap_or(false);

    if !is_valid_shape {
        set_optional_string_field(section, field_name, None);
    } else if matches!(
        section
            .get(field_name)
            .and_then(|slot| slot.as_ref())
            .and_then(|values| values.first()),
        Some(Value::String(_))
    ) {
        set_optional_string_field(section, field_name, current_value.as_deref());
    }
}

/// Collects the species names that own dedicated composition sections.
///
/// BVP reactor documents store the canonical species list in
/// `process_conditions.substances`. Those names can also appear as their own
/// sections when the document carries atomic composition data.
fn bvp_species_section_names(document: &DocumentMap) -> Vec<String> {
    let Some(process_conditions) = document.get("process_conditions") else {
        return Vec::new();
    };

    let Some(substances) = process_conditions
        .get("substances")
        .and_then(|slot| slot.as_ref())
    else {
        return Vec::new();
    };

    substances
        .iter()
        .filter_map(|value| value.as_string().cloned())
        .collect()
}

/// Renders `Option<String>` solver settings without exposing the generic optional editor.
///
/// The lower solver still accepts `linear_sys_method`, but for the BVP routes
/// used here the normal value is `None`: Banded and Sparse backends pick their
/// own linear algebra policy. If a compatibility override is needed, the GUI
/// offers only the parser-supported textual choices.
fn render_solver_backend_optional_string_choice_row(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    display_name: &str,
    choices: &[(&str, Option<&str>)],
    default_label: &str,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    normalize_optional_string_field(section, field_name);
    let current_value = optional_string_field_value(section, field_name);
    let current_label = choices
        .iter()
        .find(|(_, value)| value.map(str::to_string) == current_value)
        .map(|(label, _)| *label)
        .unwrap_or(default_label);

    let label_response = ui.label(display_name);
    if label_response.hovered() {
        if let Some(map) = help_map {
            let help_key = format!("solver_settings.{}", field_name);
            if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == help_key) {
                *current_help = help_text.to_string();
            }
        }
    }

    let mut selected_label = current_label.to_string();
    let combo_response =
        egui::ComboBox::from_id_salt(format!("solver_settings.{}.{}", display_name, field_name))
            .selected_text(selected_label.as_str())
            .show_ui(ui, |ui| {
                for (label, _) in choices {
                    ui.selectable_value(&mut selected_label, (*label).to_string(), *label);
                }
            });
    let _ = combo_response.response.labelled_by(label_response.id);

    if selected_label != current_label {
        let selected_value = choices
            .iter()
            .find(|(label, _)| *label == selected_label)
            .and_then(|(_, value)| *value);
        set_optional_string_field(section, field_name, selected_value);
        info!(
            "Field '{}' changed to: {}",
            field_name,
            selected_value.unwrap_or("None")
        );
    }

    ui.end_row();
}

/// Renders one row in the structured solver backend panel.
fn render_solver_backend_row(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    display_name: &str,
    default_value: Value,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    let values = ensure_field_slot(section, field_name, default_value);
    let label_response = ui.label(display_name);

    if label_response.hovered() {
        if let Some(map) = help_map {
            let help_key = format!("solver_settings.{}", field_name);
            if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == help_key) {
                *current_help = help_text.to_string();
            }
        }
    }

    match values {
        Some(values) => {
            if values.is_empty() {
                ui.label("None");
            } else {
                for value in values.iter_mut() {
                    render_value(ui, value, field_name);
                }
            }
        }
        None => {
            ui.label("None");
        }
    }

    ui.end_row();
}

/// Renders a solver backend row that must stay as `usize`.
///
/// The solver rejects floating-point and signed integer representations for
/// iteration counters, so this helper keeps the GUI contract aligned with the
/// parser contract.
fn render_solver_backend_usize_row(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    display_name: &str,
    default_value: usize,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    let current_value = {
        let values = ensure_field_slot(section, field_name, Value::Usize(default_value));
        values
            .as_ref()
            .and_then(|values| values.first())
            .map(try_value_as_usize)
            .unwrap_or(Ok(default_value))
    };

    let label_response = ui.label(display_name);
    if label_response.hovered() {
        if let Some(map) = help_map {
            let help_key = format!("solver_settings.{}", field_name);
            if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == help_key) {
                *current_help = help_text.to_string();
            }
        }
    }

    match current_value {
        Ok(mut current_value) => {
            let _ = normalize_usize_field(section, field_name);
            if ui
                .add(egui::DragValue::new(&mut current_value).speed(1))
                .changed()
            {
                section.insert(
                    field_name.to_string(),
                    Some(vec![Value::Usize(current_value)]),
                );
                info!("Field '{}' changed to: {}", field_name, current_value);
            }
        }
        Err(error) => {
            ui.horizontal(|ui| {
                ui.colored_label(egui::Color32::RED, error.to_string());
                if ui.button(format!("Reset to {default_value}")).clicked() {
                    section.insert(
                        field_name.to_string(),
                        Some(vec![Value::Usize(default_value)]),
                    );
                }
            });
        }
    }

    ui.end_row();
}

/// Renders a dropdown-backed string field in the solver backend panel.
fn render_solver_backend_choice_row(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    field_name: &str,
    display_name: &str,
    choices: &[&str],
    default_value: &str,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    let current_value = string_field_value(section, field_name).unwrap_or_else(|| {
        let default_value = default_value.to_string();
        set_string_field(section, field_name, &default_value);
        default_value
    });

    let label_response = ui.label(display_name);
    if label_response.hovered() {
        if let Some(map) = help_map {
            let help_key = format!("solver_settings.{}", field_name);
            if let Some((_, help_text)) = map.iter().find(|(key, _)| *key == help_key) {
                *current_help = help_text.to_string();
            }
        }
    }

    let mut selected_value = current_value.clone();
    let combo_response =
        egui::ComboBox::from_id_salt(format!("solver_settings.{}.{}", display_name, field_name))
            .selected_text(selected_value.as_str())
            .show_ui(ui, |ui| {
                for choice in choices {
                    ui.selectable_value(&mut selected_value, (*choice).to_string(), *choice);
                }
            });
    let _ = combo_response.response.labelled_by(label_response.id);

    if selected_value != current_value {
        set_string_field(section, field_name, &selected_value);
        info!("Field '{}' changed to: {}", field_name, selected_value);

        if field_name == "matrix_backend"
            || field_name == "aot_codegen_backend"
            || field_name == "aot_c_compiler"
        {
            normalize_bvp_solver_backend_section(section);
        }
    }

    ui.end_row();
}

/// Renders a user-facing execution-mode selector and keeps `generated_backend`
/// synchronized with the higher-level choice.
fn render_generated_backend_mode_row(
    ui: &mut egui::Ui,
    section: &mut HashMap<String, Option<Vec<Value>>>,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) -> bool {
    let current_backend = string_field_value(section, "generated_backend")
        .unwrap_or_else(|| "banded_lambdify".to_string());
    let current_mode = GuiExecutionBackend::from_generated_backend(&current_backend);
    let current_label = match current_mode {
        Some(GuiExecutionBackend::Aot) => "AOT".to_string(),
        Some(GuiExecutionBackend::Lambdify) => "Lambdify".to_string(),
        None => format!("Unsupported: {}", current_backend),
    };

    let label_response = ui.label("Residuals and Jacobian Backend");
    if label_response.hovered() {
        if let Some(map) = help_map {
            if let Some((_, help_text)) = map
                .iter()
                .find(|(key, _)| *key == "solver_settings.generated_backend")
            {
                *current_help = help_text.to_string();
            }
        }
    }

    let mut selected_mode = current_label.clone();
    let combo_response = egui::ComboBox::from_id_salt("solver_settings.generated_backend.mode")
        .selected_text(selected_mode.as_str())
        .show_ui(ui, |ui| {
            ui.selectable_value(&mut selected_mode, "Lambdify".to_string(), "Lambdify");
            ui.selectable_value(&mut selected_mode, "AOT".to_string(), "AOT");
        });
    let _ = combo_response.response.labelled_by(label_response.id);

    let changed = selected_mode != current_label;
    if changed {
        let matrix_backend =
            string_field_value(section, "matrix_backend").unwrap_or_else(|| "Banded".to_string());
        let backend_name = match (selected_mode.as_str(), matrix_backend.as_str()) {
            ("AOT", "Sparse") => "sparse_aot_tcc",
            ("AOT", _) => "banded_aot_tcc",
            ("Lambdify", "Sparse") => "sparse_lambdify",
            _ => "banded_lambdify",
        };
        set_string_field(section, "generated_backend", backend_name);
        normalize_bvp_solver_backend_section(section);
        info!("Field 'generated_backend' changed to: {}", backend_name);
    }

    ui.end_row();
    selected_mode == "AOT"
}

/// Renders the structured solver backend panel for BVP problems.
fn render_bvp_solver_backend_panel(
    ui: &mut egui::Ui,
    document: &mut DocumentMap,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
) {
    // Normalize persisted overrides before the widgets read the document.
    // Story tests deliberately seed partially edited solver settings, so the
    // render path must reconcile the backend contract even when the user has
    // not changed a control yet.
    ensure_bvp_solver_backend_defaults(document);
    let section = ensure_section_mut(document, "solver_settings");
    normalize_bvp_solver_backend_section(section);

    ui.separator();
    ui.label("Solver backend");
    ui.label("Default route: Lambdify + AtomView + Banded");
    ui.label("Residuals and Jacobian backend is exclusive: Lambdify or AOT.");
    ui.add_space(6.0);

    ui.label("Core solve");
    egui::Grid::new("bvp_solver_backend_core_grid")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            render_solver_backend_choice_row(
                ui,
                section,
                "scheme",
                "scheme",
                &["forward", "trapezoid"],
                "forward",
                help_map,
                current_help,
            );
            let matrix_backend_before = string_field_value(section, "matrix_backend")
                .as_deref()
                .and_then(canonical_matrix_backend)
                .unwrap_or("Banded")
                .to_string();
            render_solver_backend_choice_row(
                ui,
                section,
                "method",
                "Linear algebra backend",
                &["Banded", "Sparse"],
                "Banded",
                help_map,
                current_help,
            );
            let matrix_backend_after = string_field_value(section, "method")
                .as_deref()
                .and_then(canonical_matrix_backend)
                .unwrap_or("Banded");
            set_string_field(section, "matrix_backend", matrix_backend_after);
            if matrix_backend_after != matrix_backend_before {
                normalize_bvp_solver_backend_section(section);
            }
            render_solver_backend_choice_row(
                ui,
                section,
                "strategy",
                "strategy",
                &["Damped"],
                "Damped",
                help_map,
                current_help,
            );
            render_solver_backend_optional_string_choice_row(
                ui,
                section,
                "linear_sys_method",
                "Linear system method override",
                &[("Auto", None), ("faithful", Some("faithful"))],
                "Auto",
                help_map,
                current_help,
            );
            render_solver_backend_row(
                ui,
                section,
                "abs_tolerance",
                "abs_tolerance",
                Value::Float(1e-7),
                help_map,
                current_help,
            );
            render_solver_backend_usize_row(
                ui,
                section,
                "max_iterations",
                "max_iterations",
                100,
                help_map,
                current_help,
            );
            render_solver_backend_optional_string_choice_row(
                ui,
                section,
                "loglevel",
                "loglevel",
                &[
                    ("None", None),
                    ("error", Some("error")),
                    ("warn", Some("warn")),
                    ("info", Some("info")),
                    ("debug", Some("debug")),
                    ("trace", Some("trace")),
                ],
                "info",
                help_map,
                current_help,
            );
            render_solver_backend_row(
                ui,
                section,
                "dont_save_log",
                "dont_save_log",
                Value::Boolean(true),
                help_map,
                current_help,
            );
        });

    ui.add_space(8.0);
    ui.label("Execution and assembly");
    egui::Grid::new("bvp_solver_backend_generated_grid")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let aot_enabled =
                render_generated_backend_mode_row(ui, section, help_map, current_help);
            render_solver_backend_choice_row(
                ui,
                section,
                "symbolic_backend",
                "Symbolic backend",
                &["AtomView", "ExprLegacy"],
                "AtomView",
                help_map,
                current_help,
            );
            if aot_enabled {
                render_solver_backend_choice_row(
                    ui,
                    section,
                    "aot_codegen_backend",
                    "AOT codegen backend",
                    &["C", "Zig", "Rust"],
                    "C",
                    help_map,
                    current_help,
                );
                let codegen_is_c = string_field_value(section, "aot_codegen_backend")
                    .as_deref()
                    .and_then(GuiAotCodegenBackend::parse)
                    == Some(GuiAotCodegenBackend::C);
                if codegen_is_c {
                    render_solver_backend_choice_row(
                        ui,
                        section,
                        "aot_c_compiler",
                        "C compiler",
                        &["tcc", "gcc"],
                        "tcc",
                        help_map,
                        current_help,
                    );
                }
                render_solver_backend_choice_row(
                    ui,
                    section,
                    "aot_build_policy",
                    "AOT build policy",
                    &[
                        "build_if_missing",
                        "use_if_available",
                        "require_prebuilt",
                        "rebuild_always",
                    ],
                    "build_if_missing",
                    help_map,
                    current_help,
                );
                render_solver_backend_choice_row(
                    ui,
                    section,
                    "aot_build_profile",
                    "AOT build profile",
                    &["release", "debug"],
                    "release",
                    help_map,
                    current_help,
                );
                render_solver_backend_choice_row(
                    ui,
                    section,
                    "aot_compile_preset",
                    "AOT compile preset",
                    &["production", "fast_build", "dev_fastest"],
                    "dev_fastest",
                    help_map,
                    current_help,
                );
                render_solver_backend_choice_row(
                    ui,
                    section,
                    "aot_execution_policy",
                    "AOT execution policy",
                    &["auto", "sequential", "parallel"],
                    "auto",
                    help_map,
                    current_help,
                );
            }
            render_solver_backend_choice_row(
                ui,
                section,
                "banded_linear_solver",
                "Banded linear solver",
                &[
                    "auto",
                    "faithful",
                    "block_tridiagonal",
                    "block_tridiagonal_consistent",
                    "faer_sparse",
                    "general_partial_pivot",
                ],
                "auto",
                help_map,
                current_help,
            );
            render_solver_backend_usize_row(
                ui,
                section,
                "refinement_steps",
                "refinement_steps",
                5,
                help_map,
                current_help,
            );
            let _ = aot_enabled;
        });
}

const BVP_PHYSICS_SECTIONS: &[&str] = &[
    "process_conditions",
    "boundary_condition",
    "diffusion_coefficients",
    "reactions",
    "bounds",
];

/// Solver-level controls that are still stored in separate sections for readability.
const BVP_SOLVER_ADVANCED_SECTIONS: &[&str] = &["rel_tolerance", "strategy_params"];

const BVP_ADAPTIVE_SOLVER_SECTIONS: &[&str] = &["adaptive_strategy", "grid_refinement"];

const BVP_INITIAL_GUESS_SECTIONS: &[&str] = &["initial_guess"];

const BVP_POSTPROCESSING_SECTIONS: &[&str] = &["postprocessing"];

fn is_bvp_structured_section(section_name: &str, species_sections: &[String]) -> bool {
    section_name == "solver_settings"
        || BVP_PHYSICS_SECTIONS.contains(&section_name)
        || BVP_SOLVER_ADVANCED_SECTIONS.contains(&section_name)
        || BVP_ADAPTIVE_SOLVER_SECTIONS.contains(&section_name)
        || BVP_INITIAL_GUESS_SECTIONS.contains(&section_name)
        || BVP_POSTPROCESSING_SECTIONS.contains(&section_name)
        || species_sections
            .iter()
            .any(|species| species == section_name)
}

/// Renders a titled group of BVP sections in the structured screen layout.
///
/// The helper keeps the canonical sections visible up front while staying
/// read-only with respect to absent data: missing sections are reported, but
/// not materialized during rendering.
fn render_bvp_section_group(
    ui: &mut egui::Ui,
    title: &str,
    section_names: &[&str],
    document: &mut DocumentMap,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
    new_field_names: &mut HashMap<String, String>,
    new_field_values: &mut HashMap<String, String>,
    pending_deletion: &mut Option<DocumentDeleteRequest>,
) {
    ui.group(|ui| {
        ui.vertical(|ui| {
            ui.heading(title);
            ui.add_space(4.0);

            let mut sections_to_process = Vec::new();
            let mut sections_to_delete = Vec::new();
            let mut fields_to_delete = Vec::new();

            for section_name in section_names {
                let Some(section) = document.get_mut(*section_name) else {
                    ui.horizontal(|ui| {
                        ui.colored_label(
                            egui::Color32::from_rgb(160, 110, 20),
                            format!("Missing canonical section: {}", section_name),
                        );
                    });
                    continue;
                };

                let field_name = new_field_names
                    .entry((*section_name).to_string())
                    .or_insert_with(String::new);
                let field_value = new_field_values
                    .entry((*section_name).to_string())
                    .or_insert_with(String::new);

                let (add_field, delete_section, requested_fields) =
                    if *section_name == "grid_refinement" {
                        let (add_field, delete_section) =
                            render_grid_refinement_section(ui, section, help_map, current_help);
                        (add_field, delete_section, Vec::new())
                    } else {
                        render_section(
                            ui,
                            section_name,
                            section,
                            field_name,
                            field_value,
                            help_map,
                            current_help,
                        )
                    };

                if add_field {
                    sections_to_process.push((
                        (*section_name).to_string(),
                        field_name.clone(),
                        field_value.clone(),
                    ));
                }

                if delete_section {
                    sections_to_delete.push((*section_name).to_string());
                }
                fields_to_delete.extend(
                    requested_fields
                        .into_iter()
                        .map(|field| ((*section_name).to_string(), field)),
                );
            }

            for section_name in sections_to_delete {
                if pending_deletion.is_none() {
                    *pending_deletion = Some(DocumentDeleteRequest::Section {
                        section: section_name,
                    });
                }
            }
            for (section, field) in fields_to_delete {
                if pending_deletion.is_none() {
                    *pending_deletion = Some(DocumentDeleteRequest::Field { section, field });
                }
            }

            for (section_name, field_name, field_value) in sections_to_process {
                if !field_name.is_empty() {
                    let value = if field_value.parse::<f64>().is_ok() {
                        if let Ok(num) = field_value.parse::<f64>() {
                            Value::Float(num)
                        } else {
                            Value::String(field_value.clone())
                        }
                    } else if field_value.parse::<i64>().is_ok() {
                        if let Ok(num) = field_value.parse::<i64>() {
                            Value::Integer(num)
                        } else {
                            Value::String(field_value.clone())
                        }
                    } else if field_value.to_lowercase() == "true"
                        || field_value.to_lowercase() == "false"
                    {
                        Value::Boolean(field_value.to_lowercase() == "true")
                    } else {
                        Value::String(field_value.clone())
                    };

                    let section = document
                        .entry(section_name.clone())
                        .or_insert_with(HashMap::new);
                    match insert_field_if_absent(section, field_name.clone(), value) {
                        Ok(()) => {
                            info!("Added field '{}' to section '{}'", field_name, section_name)
                        }
                        Err(message) => warn!("{}", message),
                    }

                    new_field_names.insert(section_name.clone(), String::new());
                    new_field_values.insert(section_name, String::new());
                }
            }
        });
    });
}

/// Renders a titled group of BVP sections whose names are known only at runtime.
///
/// Species composition sections are the common example: the section names come
/// from the task document rather than a static enum-like list. Missing sections
/// are reported rather than auto-created during render.
fn render_bvp_dynamic_section_group(
    ui: &mut egui::Ui,
    title: &str,
    section_names: &[String],
    document: &mut DocumentMap,
    help_map: Option<&[(&str, &str)]>,
    current_help: &mut String,
    new_field_names: &mut HashMap<String, String>,
    new_field_values: &mut HashMap<String, String>,
    pending_deletion: &mut Option<DocumentDeleteRequest>,
) {
    if section_names.is_empty() {
        return;
    }

    ui.group(|ui| {
        ui.vertical(|ui| {
            ui.heading(title);
            ui.add_space(4.0);

            let mut sections_to_process = Vec::new();
            let mut sections_to_delete = Vec::new();
            let mut fields_to_delete = Vec::new();

            for section_name in section_names {
                let Some(section) = document.get_mut(section_name) else {
                    ui.horizontal(|ui| {
                        ui.colored_label(
                            egui::Color32::from_rgb(160, 110, 20),
                            format!("Missing canonical section: {}", section_name),
                        );
                    });
                    continue;
                };

                let field_name = new_field_names
                    .entry(section_name.clone())
                    .or_insert_with(String::new);
                let field_value = new_field_values
                    .entry(section_name.clone())
                    .or_insert_with(String::new);

                let (add_field, delete_section, requested_fields) = render_section(
                    ui,
                    section_name,
                    section,
                    field_name,
                    field_value,
                    help_map,
                    current_help,
                );

                if add_field {
                    sections_to_process.push((
                        section_name.clone(),
                        field_name.clone(),
                        field_value.clone(),
                    ));
                }

                if delete_section {
                    sections_to_delete.push(section_name.clone());
                }
                fields_to_delete.extend(
                    requested_fields
                        .into_iter()
                        .map(|field| (section_name.clone(), field)),
                );
            }

            for section_name in sections_to_delete {
                if pending_deletion.is_none() {
                    *pending_deletion = Some(DocumentDeleteRequest::Section {
                        section: section_name,
                    });
                }
            }
            for (section, field) in fields_to_delete {
                if pending_deletion.is_none() {
                    *pending_deletion = Some(DocumentDeleteRequest::Field { section, field });
                }
            }

            for (section_name, field_name, field_value) in sections_to_process {
                if !field_name.is_empty() {
                    let value = if field_value.parse::<f64>().is_ok() {
                        if let Ok(num) = field_value.parse::<f64>() {
                            Value::Float(num)
                        } else {
                            Value::String(field_value.clone())
                        }
                    } else if field_value.parse::<i64>().is_ok() {
                        if let Ok(num) = field_value.parse::<i64>() {
                            Value::Integer(num)
                        } else {
                            Value::String(field_value.clone())
                        }
                    } else if field_value.to_lowercase() == "true"
                        || field_value.to_lowercase() == "false"
                    {
                        Value::Boolean(field_value.to_lowercase() == "true")
                    } else {
                        Value::String(field_value.clone())
                    };

                    let section = document
                        .entry(section_name.clone())
                        .or_insert_with(HashMap::new);
                    match insert_field_if_absent(section, field_name.clone(), value) {
                        Ok(()) => {
                            info!("Added field '{}' to section '{}'", field_name, section_name)
                        }
                        Err(message) => warn!("{}", message),
                    }

                    new_field_names.insert(section_name.clone(), String::new());
                    new_field_values.insert(section_name, String::new());
                }
            }
        });
    });
}

/// User-visible lifecycle of one GUI calculation request.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum CalculationState {
    #[default]
    Idle,
    Running {
        run_id: u64,
    },
    Cancelling {
        run_id: u64,
    },
    Completed,
    Failed,
}

impl CalculationState {
    fn active_run_id(self) -> Option<u64> {
        match self {
            Self::Running { run_id } | Self::Cancelling { run_id } => Some(run_id),
            Self::Idle | Self::Completed | Self::Failed => None,
        }
    }

    fn is_active(self) -> bool {
        self.active_run_id().is_some()
    }
}

/// Owned plotting payload that can safely cross the worker-thread boundary.
#[derive(Debug)]
struct BvpPlotData {
    arg: String,
    values: Vec<String>,
    x_mesh: nalgebra::DVector<f64>,
    solution: nalgebra::DMatrix<f64>,
}

#[derive(Debug)]
enum BvpWorkerOutcome {
    Completed { plot: Option<BvpPlotData> },
    Failed(String),
    Cancelled,
}

#[derive(Debug)]
struct BvpWorkerMessage {
    run_id: u64,
    outcome: BvpWorkerOutcome,
}

/// Cached, user-visible result of validating the current BVP document.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct BvpGuiValidationReport {
    pub issues: Vec<String>,
}

impl BvpGuiValidationReport {
    pub(crate) fn is_valid(&self) -> bool {
        self.issues.is_empty()
    }

    fn failure_message(&self) -> String {
        format!(
            "Calculation failed: task validation failed:\n- {}",
            self.issues.join("\n- ")
        )
    }
}

/// One row in the BVP preview table.
///
/// The preview keeps the normalized document readable without mutating the
/// loaded task state or inventing a new data model.
#[derive(Clone, Debug, PartialEq, Eq, Tabled)]
pub(crate) struct BvpPreviewParameterRow {
    #[tabled(rename = "Section")]
    pub(crate) section: String,
    #[tabled(rename = "Field")]
    pub(crate) field: String,
    #[tabled(rename = "Value")]
    pub(crate) value: String,
}

/// One row in the equation preview table.
#[derive(Clone, Debug, PartialEq, Eq, Tabled)]
pub(crate) struct BvpPreviewEquationRow {
    #[tabled(rename = "#")]
    pub(crate) index: usize,
    #[tabled(rename = "Unknown")]
    pub(crate) unknown: String,
    #[tabled(rename = "Equation")]
    pub(crate) equation: String,
}

/// Read-only snapshot used by the BVP preview action.
///
/// The snapshot is intentionally decoupled from egui so tests can inspect the
/// table rows directly without capturing terminal output.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct BvpTaskPreviewSnapshot {
    pub(crate) validation_report: BvpGuiValidationReport,
    pub(crate) parameter_rows: Vec<BvpPreviewParameterRow>,
    pub(crate) equation_rows: Vec<BvpPreviewEquationRow>,
    pub(crate) equation_build_error: Option<String>,
}

/// Validate and canonicalize one GUI document without starting solver work.
fn validate_bvp_gui_document(
    document: &DocumentMap,
    problem: &ProblemsEnum,
) -> (DocumentMap, BvpGuiValidationReport) {
    let mut normalized = document.clone();
    let mut report = BvpGuiValidationReport::default();

    if let Err(error) = normalize_bvp_usize_fields(&mut normalized) {
        report.issues.push(error);
    }
    if let Err(error) = validate_aot_toolchain_fields(&normalized) {
        report.issues.push(format!("AOT security: {error}"));
    }

    match problem {
        ProblemsEnum::BVPSimple => {
            let reactor_report = validate_reactor_bvp_task_map(&normalized);
            report.issues.extend(
                reactor_report
                    .issues
                    .into_iter()
                    .map(|issue| format!("{}: {}", issue.stage, issue.message)),
            );
        }
        ProblemsEnum::None => report
            .issues
            .push("No calculation is available for the current problem type".to_string()),
    }

    (normalized, report)
}

fn bvp_preview_section_rank(section_name: &str, species_sections: &[String]) -> usize {
    if section_name == "solver_settings" {
        0
    } else if BVP_PHYSICS_SECTIONS.contains(&section_name) {
        1
    } else if species_sections
        .iter()
        .any(|species| species == section_name)
    {
        2
    } else if BVP_SOLVER_ADVANCED_SECTIONS.contains(&section_name) {
        3
    } else if BVP_ADAPTIVE_SOLVER_SECTIONS.contains(&section_name) {
        4
    } else if BVP_INITIAL_GUESS_SECTIONS.contains(&section_name) {
        5
    } else if BVP_POSTPROCESSING_SECTIONS.contains(&section_name) {
        6
    } else {
        7
    }
}

/// Collects the normalized BVP document into a table-friendly snapshot.
///
/// The preview keeps raw compatibility sections visible but sorts the canonical
/// BVP groups first so users can scan the document in a predictable order.
pub(crate) fn collect_bvp_preview_parameter_rows(
    document: &DocumentMap,
) -> Vec<BvpPreviewParameterRow> {
    let species_sections = bvp_species_section_names(document);
    let mut section_names = document.keys().cloned().collect::<Vec<_>>();
    section_names.sort_by_key(|section_name| {
        (
            bvp_preview_section_rank(section_name, &species_sections),
            section_name.clone(),
        )
    });

    let mut rows = Vec::new();
    for section_name in section_names {
        let Some(section) = document.get(&section_name) else {
            continue;
        };
        let mut field_names = section.keys().cloned().collect::<Vec<_>>();
        field_names.sort();
        for field_name in field_names {
            let value = match section.get(&field_name).and_then(|slot| slot.as_ref()) {
                Some(values) if !values.is_empty() => values_to_document_string(values),
                Some(_) => "<empty>".to_string(),
                None => "None".to_string(),
            };
            rows.push(BvpPreviewParameterRow {
                section: section_name.clone(),
                field: field_name,
                value,
            });
        }
    }
    rows
}

/// Collects the assembled BVP equations into a stable table snapshot.
///
/// The row order follows the solver order so the preview matches the numerical
/// system that will be handed to RustedSciThe.
pub(crate) fn collect_bvp_preview_equation_rows(
    reactor: &SimpleReactorTask,
) -> Vec<BvpPreviewEquationRow> {
    let equation_count = reactor
        .solver
        .eq_system
        .len()
        .min(reactor.solver.unknowns.len());

    let mut rows = Vec::with_capacity(equation_count);
    for index in 0..equation_count {
        let unknown = reactor
            .solver
            .unknowns
            .get(index)
            .cloned()
            .unwrap_or_else(|| format!("eq_{}", index + 1));
        let equation = reactor.solver.eq_system[index].pretty_print();
        rows.push(BvpPreviewEquationRow {
            index: index + 1,
            unknown,
            equation,
        });
    }
    rows
}

/// Builds a read-only preview snapshot for the current BVP task.
///
/// The builder uses the same validation and equation assembly path as a normal
/// solve request, but it runs on a local clone and never mutates GUI state.
pub(crate) fn build_bvp_task_preview_snapshot(
    document: &DocumentMap,
    problem: &ProblemsEnum,
) -> BvpTaskPreviewSnapshot {
    let (normalized_document, validation_report) = validate_bvp_gui_document(document, problem);
    let parameter_rows = collect_bvp_preview_parameter_rows(&normalized_document);
    let mut snapshot = BvpTaskPreviewSnapshot {
        validation_report,
        parameter_rows,
        equation_rows: Vec::new(),
        equation_build_error: None,
    };

    match problem {
        ProblemsEnum::BVPSimple => {
            let mut reactor = SimpleReactorTask::new();
            if let Err(error) = reactor.set_reactor_params_from_hashmap(&normalized_document) {
                snapshot.equation_build_error =
                    Some(format!("Preview could not prepare the BVP task: {}", error));
            } else if let Err(error) = reactor.setup_bvp() {
                snapshot.equation_build_error = Some(format!(
                    "Preview could not assemble BVP equations: {}",
                    error
                ));
            } else {
                snapshot.equation_rows = collect_bvp_preview_equation_rows(&reactor);
            }
        }
        ProblemsEnum::None => {
            snapshot.equation_build_error =
                Some("Preview is only available for the BVP problem type.".to_string());
        }
    }

    snapshot
}

/// Prints a BVP preview snapshot as a readable console report.
///
/// The output is intentionally side-effect free with respect to task data.
pub(crate) fn print_bvp_task_preview(snapshot: &BvpTaskPreviewSnapshot) {
    info!("=== BVP TASK PREVIEW ===");
    if snapshot.validation_report.is_valid() {
        info!("Validation: ok");
    } else {
        warn!("Validation issues:");
        for issue in &snapshot.validation_report.issues {
            warn!("- {}", issue);
        }
    }

    info!("--- Normalized document ---");
    if snapshot.parameter_rows.is_empty() {
        info!("(no document fields to preview)");
    } else {
        info!(
            "{}",
            Table::new(snapshot.parameter_rows.clone()).to_string()
        );
    }

    if let Some(error) = &snapshot.equation_build_error {
        warn!("{}", error);
    }

    info!("--- Assembled equations ---");
    if snapshot.equation_rows.is_empty() {
        info!("(no equations to preview)");
    } else {
        info!("{}", Table::new(snapshot.equation_rows.clone()).to_string());
    }
}

/// Execute one reactor task without touching egui state.
fn execute_bvp_calculation(mut document: DocumentMap, problem: ProblemsEnum) -> BvpWorkerOutcome {
    if let Err(message) = normalize_bvp_usize_fields(&mut document) {
        return BvpWorkerOutcome::Failed(format!("Calculation blocked: {message}"));
    }
    if let Err(message) = validate_aot_toolchain_fields(&document) {
        return BvpWorkerOutcome::Failed(format!("Calculation blocked: {}", message));
    }

    match problem {
        ProblemsEnum::BVPSimple => {
            info!("Starting BVP Simple calculation in worker thread.");
            let gui_plot = document
                .get("postprocessing")
                .and_then(|postproc| postproc.get("gui_plot"))
                .and_then(|slot| slot.as_ref())
                .and_then(|values| values.first())
                .and_then(Value::as_boolean)
                .unwrap_or(true);

            let mut reactor = SimpleReactorTask::new();
            if let Err(error) = reactor.solve_from_parsed(document) {
                return BvpWorkerOutcome::Failed(format!("Calculation failed: {}", error));
            }
            if !gui_plot {
                return BvpWorkerOutcome::Completed { plot: None };
            }

            let Some(solution) = reactor.solver.solution.take() else {
                return BvpWorkerOutcome::Failed(
                    "Calculation completed, but no solution matrix was returned.".to_string(),
                );
            };
            let Some(x_mesh) = reactor.solver.x_mesh.take() else {
                return BvpWorkerOutcome::Failed(
                    "Calculation completed, but no mesh was returned.".to_string(),
                );
            };
            let values = std::mem::take(&mut reactor.solver.unknowns);
            if x_mesh.len() != solution.nrows() {
                return BvpWorkerOutcome::Failed(format!(
                    "Solver result shape mismatch: mesh has {} points but solution has {} rows",
                    x_mesh.len(),
                    solution.nrows()
                ));
            }
            if values.len() != solution.ncols() {
                return BvpWorkerOutcome::Failed(format!(
                    "Solver result shape mismatch: {} unknowns but solution has {} columns",
                    values.len(),
                    solution.ncols()
                ));
            }

            BvpWorkerOutcome::Completed {
                plot: Some(BvpPlotData {
                    arg: reactor.solver.arg_name.clone(),
                    values,
                    x_mesh,
                    solution,
                }),
            }
        }
        ProblemsEnum::None => BvpWorkerOutcome::Failed(
            "No calculation is available for the current problem type.".to_string(),
        ),
    }
}

/// Main application state for the combustion reactor GUI.
///
/// This struct maintains all the state necessary for the GUI application,
/// including the current document being edited, UI state for adding new
/// fields and sections, file management, and window visibility flags.
///
/// # Fields
///
/// * `document` - The parsed configuration document as a nested HashMap structure
/// * `new_section_name` - Input buffer for creating new sections
/// * `new_field_names` - Per-section input buffers for new field names
/// * `new_field_values` - Per-section input buffers for new field values
/// * `selected_problem` - Currently selected problem type
/// * `current_file_path` - Path to the currently loaded/saved file (if any)
/// * `document_lifecycle` - Typed file lifecycle state and latest load/save result
/// * `show_help_window` - Flag controlling help window visibility
///
/// # Design Notes
///
/// The application uses a document-centric approach where all configuration
/// data is stored in a DocumentMap structure that can be serialized to/from
/// text files and passed directly to solvers.
pub struct CombustionApp {
    pub document: DocumentMap,
    pub new_section_name: String,
    pub new_field_names: HashMap<String, String>,
    pub new_field_values: HashMap<String, String>,
    pub selected_problem: ProblemsEnum,
    pub current_file_path: Option<std::path::PathBuf>,
    pub(crate) document_lifecycle: DocumentLifecycleState,
    pub show_help_window: bool,
    pub current_help_info: String,
    pub last_run_message: Option<String>,
    pub last_run_is_error: bool,
    pub plot_window: Option<PlotWindow>,
    pub(crate) pending_aot_confirmation: Option<AotToolchainRequest>,
    pub(crate) pending_document_load: Option<std::path::PathBuf>,
    pub(crate) pending_document_deletion: Option<DocumentDeleteRequest>,
    pub(crate) pending_window_close: bool,
    pub(crate) pending_problem_switch: Option<ProblemsEnum>,
    confirmed_aot_request: Option<AotToolchainRequest>,
    pub(crate) calculation_state: CalculationState,
    calculation_receiver: Option<mpsc::Receiver<BvpWorkerMessage>>,
    calculation_cancel: Option<Arc<AtomicBool>>,
    next_run_id: u64,
    validation_fingerprint: Option<String>,
    pub(crate) validation_report: BvpGuiValidationReport,
    pub(crate) bvp_gui_config: BvpGuiConfig,
    pub(crate) bvp_gui_migration_report: BvpGuiMigrationReport,
}

impl CombustionApp {
    /// Creates a new CombustionApp instance with default problem type.
    ///
    /// This is a convenience constructor that creates an app with the default
    /// problem type (None), resulting in an empty document ready for user input.
    pub fn new() -> Self {
        Self::new_with_problem(ProblemsEnum::default())
    }

    /// Converts the current document to a string representation suitable for file storage.
    ///
    /// This method serializes the DocumentMap structure back into the text format
    /// that can be parsed by DocumentParser. It handles all Value types and
    /// maintains the section-field hierarchy.
    ///
    /// # Returns
    ///
    /// A string representation of the document in the format:
    /// ```text
    /// section_name
    /// field_name: value
    /// field_name2: value2
    /// ```
    ///
    /// # Value Type Handling
    ///
    /// * Primitive types (String, Float, Integer, etc.) are formatted directly
    /// * Vectors are formatted as comma-separated lists in brackets
    /// * Optional values show "Some(value)" or "None"
    /// * None fields show "field_name: None"
    pub fn document_to_string(&self) -> String {
        let mut result = String::new();
        let mut section_names = self.document.keys().collect::<Vec<_>>();
        section_names.sort();

        for section_name in section_names {
            let section_data = &self.document[section_name];
            if section_data.is_empty() {
                continue;
            }
            result.push_str(section_name);
            result.push('\n');

            let mut field_names = section_data.keys().collect::<Vec<_>>();
            field_names.sort();

            for field_name in field_names {
                let maybe_values = &section_data[field_name];
                if let Some(values) = maybe_values {
                    result.push_str(&format!(
                        "{}: {}\n",
                        field_name,
                        values_to_document_string(values)
                    ));
                } else {
                    result.push_str(&format!("{}: None\n", field_name));
                }
            }
        }
        result
    }

    /// Renders the typed BVP editor against a normalized working snapshot.
    ///
    /// The visible editor works from the typed configuration, then commits the
    /// edited snapshot back to the live document once the structured pass ends.
    fn render_bvp_structured_editor(
        &mut self,
        ui: &mut egui::Ui,
        help_map: Option<&[(&str, &str)]>,
    ) {
        let mut editor_document = self.document.clone();
        render_bvp_solver_backend_panel(
            ui,
            &mut editor_document,
            help_map,
            &mut self.current_help_info,
        );
        render_bvp_section_group(
            ui,
            "Advanced solver settings",
            BVP_SOLVER_ADVANCED_SECTIONS,
            &mut editor_document,
            help_map,
            &mut self.current_help_info,
            &mut self.new_field_names,
            &mut self.new_field_values,
            &mut self.pending_document_deletion,
        );
        render_bvp_adaptive_grid_panel(
            ui,
            &mut editor_document,
            help_map,
            &mut self.current_help_info,
        );
        render_bvp_section_group(
            ui,
            "Initial guess",
            BVP_INITIAL_GUESS_SECTIONS,
            &mut editor_document,
            help_map,
            &mut self.current_help_info,
            &mut self.new_field_names,
            &mut self.new_field_values,
            &mut self.pending_document_deletion,
        );
        render_bvp_section_group(
            ui,
            "Postprocessing",
            BVP_POSTPROCESSING_SECTIONS,
            &mut editor_document,
            help_map,
            &mut self.current_help_info,
            &mut self.new_field_names,
            &mut self.new_field_values,
            &mut self.pending_document_deletion,
        );

        self.document = editor_document;
        self.sync_document_state_after_edit();
    }

    /// Render confirmation for deleting structured BVP data.
    fn show_document_deletion_confirmation(&mut self, ctx: &egui::Context) {
        let Some(request) = self.pending_document_deletion.clone() else {
            return;
        };

        let mut confirm = false;
        let mut cancel = false;
        egui::Window::new("Confirm BVP data deletion")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(format!(
                    "Delete {}? This may make the reactor task incomplete.",
                    request.description()
                ));
                ui.horizontal(|ui| {
                    if ui.button("Cancel deletion").clicked() {
                        cancel = true;
                    }
                    if ui.button("Confirm deletion").clicked() {
                        confirm = true;
                    }
                });
            });

        if cancel {
            self.pending_document_deletion = None;
            self.last_run_message =
                Some("Deletion cancelled; document data was preserved.".to_string());
            self.last_run_is_error = false;
        } else if confirm {
            let removed = apply_document_deletion(&mut self.document, &request);
            self.pending_document_deletion = None;
            if removed {
                self.last_run_message = Some(format!("Deleted {}.", request.description()));
                self.last_run_is_error = false;
            } else {
                self.last_run_message = Some(format!(
                    "Could not delete {}; it no longer exists.",
                    request.description()
                ));
                self.last_run_is_error = true;
            }
            self.sync_document_state_after_edit();
        }
    }

    /// Render confirmation for replacing the current task via Read Task.
    fn show_document_load_confirmation(&mut self, ctx: &egui::Context) {
        let Some(path) = self.pending_document_load.clone() else {
            return;
        };

        let mut confirm = false;
        let mut cancel = false;
        egui::Window::new("Confirm task load")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(format!(
                    "Load {}? This will replace the current task state.",
                    path.display()
                ));
                ui.label("Confirm only if you want to discard unsaved edits or running work.");
                ui.horizontal(|ui| {
                    if ui.button("Cancel load").clicked() {
                        cancel = true;
                    }
                    if ui.button("Load task").clicked() {
                        confirm = true;
                    }
                });
            });

        if cancel {
            self.cancel_pending_document_load();
        } else if confirm {
            let _ = self.confirm_pending_document_load();
        }
    }

    /// Render confirmation for closing the main window with unsaved work.
    fn show_window_close_confirmation(&mut self, ctx: &egui::Context, open: &mut bool) {
        if !self.pending_window_close {
            return;
        }

        let mut confirm = false;
        let mut cancel = false;
        egui::Window::new("Confirm window close")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label("Close the window and discard the current task state?");
                ui.label("Confirm only if you do not need the unsaved changes anymore.");
                ui.horizontal(|ui| {
                    if ui.button("Cancel close").clicked() {
                        cancel = true;
                    }
                    if ui.button("Close window").clicked() {
                        confirm = true;
                    }
                });
            });

        if cancel {
            self.cancel_pending_window_close();
            *open = true;
        } else if confirm {
            if self.confirm_pending_window_close() {
                *open = false;
            }
        } else {
            *open = true;
        }
    }

    /// Render confirmation for switching the active problem type.
    fn show_problem_switch_confirmation(&mut self, ctx: &egui::Context) {
        let Some(requested_problem) = self.pending_problem_switch.clone() else {
            return;
        };

        let mut confirm = false;
        let mut cancel = false;
        egui::Window::new("Confirm problem switch")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(format!(
                    "Switching to {} will replace the current task state.",
                    requested_problem
                ));
                ui.label("Confirm only if you want to discard the current edits or running work.");
                ui.horizontal(|ui| {
                    if ui.button("Cancel switch").clicked() {
                        cancel = true;
                    }
                    if ui.button("Switch problem").clicked() {
                        confirm = true;
                    }
                });
            });

        if cancel {
            self.cancel_pending_problem_switch();
        } else if confirm {
            if self.confirm_pending_problem_switch() {
                info!("Switched to problem type: {}", self.selected_problem);
            }
        }
    }

    /// Creates a new CombustionApp instance with a specific problem type.
    ///
    /// This constructor initializes the application with a template document
    /// corresponding to the specified problem type. If a template is available,
    /// it is parsed and loaded into the document field.
    ///
    /// # Arguments
    ///
    /// * `problem` - The problem type to initialize with
    ///
    /// # Template Loading
    ///
    /// If the problem type has an associated template:
    /// 1. The template string is retrieved
    /// 2. A DocumentParser is created and used to parse the template
    /// 3. The parsed result is stored in the document field
    /// 4. Parse errors are logged but don't prevent app creation
    ///
    /// If no template is available, an empty document is created.
    pub fn new_with_problem(problem: ProblemsEnum) -> Self {
        let document = match problem.get_problem() {
            Some(template) => {
                let mut parser = DocumentParser::new(template.to_string());
                match parser.parse_document() {
                    Ok(_) => info!("Parsed BVP template successfully"),
                    Err(e) => warn!("Error parsing BVP template: {}", e),
                }
                parser.get_result().cloned().unwrap_or_else(HashMap::new)
            }
            None => HashMap::new(),
        };
        let (bvp_gui_config, bvp_gui_migration_report) = BvpGuiConfig::from_document(&document);

        Self {
            document,
            new_section_name: String::new(),
            new_field_names: HashMap::new(),
            new_field_values: HashMap::new(),
            selected_problem: problem,
            current_file_path: None,
            document_lifecycle: DocumentLifecycleState::default(),
            show_help_window: false,
            current_help_info: String::new(),
            last_run_message: None,
            last_run_is_error: false,
            plot_window: None,
            pending_aot_confirmation: None,
            pending_document_load: None,
            pending_document_deletion: None,
            pending_window_close: false,
            pending_problem_switch: None,
            confirmed_aot_request: None,
            calculation_state: CalculationState::Idle,
            calculation_receiver: None,
            calculation_cancel: None,
            next_run_id: 1,
            validation_fingerprint: None,
            validation_report: BvpGuiValidationReport::default(),
            bvp_gui_config,
            bvp_gui_migration_report,
        }
        .with_seeded_defaults()
    }

    /// Seeds problem-specific defaults after constructing the app state.
    fn with_seeded_defaults(mut self) -> Self {
        if matches!(self.selected_problem, ProblemsEnum::BVPSimple) {
            seed_bvp_document_defaults(&mut self.document);
        }
        self.document_lifecycle
            .mark_clean_snapshot(self.current_validation_fingerprint());
        self.refresh_bvp_gui_snapshot();
        self
    }

    /// Updates the active file path in both the compatibility field and the
    /// explicit lifecycle state.
    /// Renders the main application window and handles all user interactions.
    ///
    /// This is the primary UI method that creates and manages the main application
    /// window, including all menus, controls, and content areas. It handles:
    ///
    /// * File operations (Save, Save As, Read Task)
    /// * Problem type selection and template switching
    /// * Document section and field editing
    /// * Help window management
    /// * Calculation execution
    ///
    /// # Arguments
    ///
    /// * `ctx` - The egui context for rendering
    /// * `open` - Mutable reference to window open state
    ///
    /// # UI Layout
    ///
    /// The window is organized into several sections:
    /// 1. Top toolbar with File menu, problem selector, and Read Task button
    /// 2. Main content area with scrollable document sections
    /// 3. Section creation controls
    /// 4. Bottom calculation button
    /// 5. Conditional help window overlay
    ///
    /// # State Management
    ///
    /// The method handles complex state updates including:
    /// * Problem type changes (reloads document template)
    /// * File path tracking for save operations
    /// * Dynamic field and section addition/deletion
    /// * Input validation and type conversion
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        self.poll_calculation_worker(ctx);
        self.refresh_validation_report();
        self.refresh_bvp_gui_snapshot();
        let was_open = *open;
        let mut window_open = *open;
        egui::Window::new("Gas-phase Combustion/Steady State Plug Flow")
            .open(&mut window_open)
            .default_size(ctx.content_rect().size())
            .resizable(true)
            .collapsible(true)
            .show(ctx, |ui| {
                // Top bar with file menu, dropdown and read task button
                ui.horizontal(|ui| {
                    // File menu dropdown on the far left
                    egui::ComboBox::from_label("File")
                        .selected_text("File")
                        .show_ui(ui, |ui| {
                            if ui.button("💾 Save").clicked() {
                                if let Some(path) = &self.current_file_path {
                                    let _ = self.save_document(path.clone());
                                } else {
                                    if let Some(path) = rfd::FileDialog::new()
                                        .add_filter("text", &["txt"])
                                        .save_file()
                                    {
                                        let _ = self.save_document(path);
                                    }
                                }
                            }
                            if ui.button("💾 Save As...").clicked() {
                                if let Some(path) = rfd::FileDialog::new()
                                    .add_filter("text", &["txt"])
                                    .save_file()
                                {
                                    let _ = self.save_document(path);
                                }
                            }
                            if ui.button("❓ Help").clicked() {
                                self.show_help_window = true;
                            }
                        });

                    ui.separator();

                    // Problem selection dropdown
                    ui.label("Problem Type:");
                    let current_problem = self.selected_problem.clone();
                    let mut ui_selected_problem = current_problem.clone();
                    egui::ComboBox::from_label("")
                        .selected_text(format!("{}", current_problem))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut ui_selected_problem, ProblemsEnum::None, "None");
                            ui.selectable_value(
                                &mut ui_selected_problem,
                                ProblemsEnum::BVPSimple,
                                "BVP Simple",
                            );
                        });

                    // If problem changed, either switch immediately or queue confirmation.
                    if ui_selected_problem != current_problem {
                        let _ = self.request_problem_switch(ui_selected_problem);
                    }

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("📁 Read Task").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("text", &["txt"])
                                .add_filter("all", &["*"])
                                .pick_file()
                            {
                                let _ = self.request_document_load(path);
                            }
                        }
                });
            });

        ui.separator();
        ui.add_space(10.0);

                // Main content area
                ui.label("Gas-phase Combustion and Steady State Plug Flow Analysis");

                // Context-sensitive help display (always present to prevent layout shifts)
                let help_map = match self.selected_problem {
                    ProblemsEnum::BVPSimple => Some(FIELD_HELP_MAP_BVP),
                    ProblemsEnum::None => None,
                };

                if matches!(self.selected_problem, ProblemsEnum::BVPSimple) {
                    self.render_bvp_structured_editor(ui, help_map);
                }

                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("ℹ");
                    ui.style_mut().text_styles.insert(
                        egui::TextStyle::Body,
                        egui::FontId::new(14.0, egui::FontFamily::Proportional),
                    );
                    let help_width = (ui.available_width() - 30.0).max(0.0);
                    let display_text = if self.current_help_info.is_empty() {
                        "Hover over fields and sections for help information"
                    } else {
                        &self.current_help_info
                    };
                    ui.add_sized(
                        [help_width, 28.0],
                        egui::Label::new(
                            egui::RichText::new(display_text)
                                .color(egui::Color32::BLACK)
                                .size(14.0),
                        ),
                    );
                });
                if !self.bvp_gui_migration_report.warnings.is_empty() {
                    ui.group(|ui| {
                        ui.colored_label(
                            egui::Color32::from_rgb(150, 110, 20),
                            "Compatibility notes",
                        );
                        for warning in &self.bvp_gui_migration_report.warnings {
                            ui.label(format!("- {warning}"));
                        }
                    });
                    ui.add_space(8.0);
                }
                ui.add_space(20.0);

                // Render document sections
                egui::ScrollArea::vertical().show(ui, |ui| {
                    if matches!(self.selected_problem, ProblemsEnum::BVPSimple) {
                        let bvp_species_sections = bvp_species_section_names(&self.document);
                        render_bvp_section_group(
                            ui,
                            "Physics",
                            BVP_PHYSICS_SECTIONS,
                            &mut self.document,
                            help_map,
                            &mut self.current_help_info,
                            &mut self.new_field_names,
                            &mut self.new_field_values,
                            &mut self.pending_document_deletion,
                        );
                        ui.add_space(8.0);

                        render_bvp_dynamic_section_group(
                            ui,
                            "Atomic composition",
                            &bvp_species_sections,
                            &mut self.document,
                            help_map,
                            &mut self.current_help_info,
                            &mut self.new_field_names,
                            &mut self.new_field_values,
                            &mut self.pending_document_deletion,
                        );
                        ui.add_space(8.0);

                        ui.add_space(12.0);
                        egui::CollapsingHeader::new("Legacy raw document")
                            .default_open(false)
                            .show(ui, |ui| {
                                let mut sections_to_process = Vec::new();
                                let mut sections_to_delete = Vec::new();
                                let mut fields_to_delete = Vec::new();

                                let section_names = sorted_document_section_names(&self.document);
                                for section_name in section_names {
                                    if is_bvp_structured_section(section_name.as_str(), &bvp_species_sections)
                                    {
                                        continue;
                                    }

                                    let Some(section_data) = self.document.get_mut(&section_name) else {
                                        continue;
                                    };

                                    let field_name = self
                                        .new_field_names
                                        .entry(section_name.clone())
                                        .or_insert_with(String::new);
                                    let field_value = self
                                        .new_field_values
                                        .entry(section_name.clone())
                                        .or_insert_with(String::new);

                                    let (add_field, delete_section, requested_fields) = render_section(
                                        ui,
                                        &section_name,
                                        section_data,
                                        field_name,
                                        field_value,
                                        help_map,
                                        &mut self.current_help_info,
                                    );

                                    if add_field {
                                        sections_to_process.push((
                                            section_name.clone(),
                                            field_name.clone(),
                                            field_value.clone(),
                                        ));
                                    }

                                    if delete_section {
                                        sections_to_delete.push(section_name.clone());
                                    }
                                    fields_to_delete.extend(
                                        requested_fields
                                            .into_iter()
                                            .map(|field| (section_name.clone(), field)),
                                    );
                                }

                                for section_name in sections_to_delete {
                                    self.document.remove(&section_name);
                                    info!("Deleted section: {}", section_name);
                                }
                                for (section, field) in fields_to_delete {
                                    let request = DocumentDeleteRequest::Field { section, field };
                                    if apply_document_deletion(&mut self.document, &request) {
                                        info!("Deleted optional {}", request.description());
                                    }
                                }

                                for (section_name, field_name, field_value) in sections_to_process {
                                    if !field_name.is_empty() {
                                        let value = if field_value.parse::<f64>().is_ok() {
                                            Value::Float(field_value.parse().unwrap())
                                        } else if field_value.parse::<i64>().is_ok() {
                                            Value::Integer(field_value.parse().unwrap())
                                        } else if field_value == "true"
                                            || field_value == "false"
                                        {
                                            Value::Boolean(field_value.parse().unwrap())
                                        } else {
                                            Value::String(field_value.clone())
                                        };

                                        if let Some(section) = self.document.get_mut(&section_name) {
                                            match insert_field_if_absent(
                                                section,
                                                field_name.clone(),
                                                value,
                                            ) {
                                                Ok(()) => info!(
                                                    "Added field '{}' to section '{}' with value: {}",
                                                    field_name, section_name, field_value
                                                ),
                                                Err(message) => warn!("{}", message),
                                            }
                                        }

                                        self.new_field_names
                                            .insert(section_name.clone(), String::new());
                                        self.new_field_values
                                            .insert(section_name, String::new());
                                    }
                                }
                            });
                    } else {
                        let mut sections_to_process = Vec::new();
                        let mut sections_to_delete = Vec::new();
                        let mut fields_to_delete = Vec::new();

                        let section_names = sorted_document_section_names(&self.document);
                        for section_name in section_names {
                            let Some(section_data) = self.document.get_mut(&section_name) else {
                                continue;
                            };

                            let field_name = self
                                .new_field_names
                                .entry(section_name.clone())
                                .or_insert_with(String::new);
                            let field_value = self
                                .new_field_values
                                .entry(section_name.clone())
                                .or_insert_with(String::new);

                            let (add_field, delete_section, requested_fields) = render_section(
                                ui,
                                &section_name,
                                section_data,
                                field_name,
                                field_value,
                                help_map,
                                &mut self.current_help_info,
                            );

                            if add_field {
                                sections_to_process.push((
                                    section_name.clone(),
                                    field_name.clone(),
                                    field_value.clone(),
                                ));
                            }

                            if delete_section {
                                sections_to_delete.push(section_name.clone());
                            }
                            fields_to_delete.extend(
                                requested_fields
                                    .into_iter()
                                    .map(|field| (section_name.clone(), field)),
                            );
                        }

                        // Delete sections outside the iterator
                        for section_name in sections_to_delete {
                            self.document.remove(&section_name);
                            info!("Deleted section: {}", section_name);
                        }
                        for (section, field) in fields_to_delete {
                            let request = DocumentDeleteRequest::Field { section, field };
                            if apply_document_deletion(&mut self.document, &request) {
                                info!("Deleted optional {}", request.description());
                            }
                        }

                        // Process field additions
                        for (section_name, field_name, field_value) in sections_to_process {
                            if !field_name.is_empty() {
                                let value = if field_value.parse::<f64>().is_ok() {
                                    Value::Float(field_value.parse().unwrap())
                                } else if field_value.parse::<i64>().is_ok() {
                                    Value::Integer(field_value.parse().unwrap())
                                } else if field_value == "true" || field_value == "false" {
                                    Value::Boolean(field_value.parse().unwrap())
                                } else {
                                    Value::String(field_value.clone())
                                };

                                if let Some(section) = self.document.get_mut(&section_name) {
                                    match insert_field_if_absent(
                                        section,
                                        field_name.clone(),
                                        value,
                                    ) {
                                        Ok(()) => info!(
                                            "Added field '{}' to section '{}' with value: {}",
                                            field_name, section_name, field_value
                                        ),
                                        Err(message) => warn!("{}", message),
                                    }
                                }

                                self.new_field_names
                                    .insert(section_name.clone(), String::new());
                                self.new_field_values.insert(section_name, String::new());
                            }
                        }
                    }
                    self.sync_document_state_after_edit();

                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label("New section name:");
                        ui.text_edit_singleline(&mut self.new_section_name);
                        if ui.button("Create Section").clicked()
                            && !self.new_section_name.is_empty()
                        {
                            match insert_section_if_absent(
                                &mut self.document,
                                self.new_section_name.clone(),
                            ) {
                                Ok(()) => info!("Created new section: {}", self.new_section_name),
                                Err(message) => warn!("{}", message),
                            }
                            self.new_section_name.clear();
                            self.sync_document_state_after_edit();
                        }
                    });
                });

                // Big RUN CALCULATION button at the bottom
                ui.separator();
                ui.add_space(10.0);
                self.refresh_document_lifecycle_dirty_state();
                if self.document_lifecycle.dirty {
                    ui.colored_label(
                        egui::Color32::from_rgb(180, 120, 20),
                        "Unsaved changes",
                    );
                    ui.add_space(4.0);
                }
                if let Some(message) = self.document_lifecycle.last_status_text() {
                    let status_color = if self.document_lifecycle.last_status_is_error() {
                        egui::Color32::from_rgb(180, 30, 30)
                    } else {
                        egui::Color32::from_rgb(20, 120, 40)
                    };
                    ui.label(
                        egui::RichText::new(message)
                            .color(status_color)
                            .size(14.0),
                    );
                    ui.add_space(4.0);
                }
                if let Some(message) = &self.last_run_message {
                    let status_color = if self.last_run_is_error {
                        egui::Color32::from_rgb(180, 30, 30)
                    } else {
                        egui::Color32::from_rgb(20, 120, 40)
                    };
                    ui.label(
                        egui::RichText::new(message.as_str())
                            .color(status_color)
                            .size(14.0),
                    );
                    ui.add_space(4.0);
                }
                ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                    let calculation_active = self.calculation_state.is_active();
                    let task_valid = self.validation_report.is_valid();
                    if !task_valid {
                        ui.group(|ui| {
                            ui.colored_label(egui::Color32::from_rgb(180, 30, 30), "Task validation");
                            for issue in &self.validation_report.issues {
                                ui.label(format!("- {issue}"));
                            }
                        });
                        ui.add_space(4.0);
                    }
                    ui.horizontal_centered(|ui| {
                        if matches!(self.selected_problem, ProblemsEnum::BVPSimple) {
                            let preview_button = egui::Button::new("Preview Task")
                                .min_size(egui::Vec2::new(160.0, 50.0))
                                .fill(egui::Color32::from_rgb(140, 150, 180));
                            if ui.add_sized([160.0, 50.0], preview_button).clicked() {
                                self.preview_bvp_task();
                            }
                            ui.add_space(8.0);
                        }

                        let button = egui::Button::new("🚀 RUN CALCULATION!")
                            .min_size(egui::Vec2::new(200.0, 50.0))
                            .fill(egui::Color32::from_rgb(160, 160, 160));

                        if ui
                            .add_enabled_ui(!calculation_active && task_valid, |ui| {
                                ui.add_sized([200.0, 50.0], button)
                            })
                            .inner
                            .clicked()
                        {
                            self.request_calculation();
                        }
                    });
                    if calculation_active {
                        let cancelling = matches!(
                            self.calculation_state,
                            CalculationState::Cancelling { .. }
                        );
                        if ui
                            .add_enabled(!cancelling, egui::Button::new("Cancel calculation"))
                            .clicked()
                        {
                            self.cancel_calculation();
                        }
                    }
                });
            });

        if was_open && !window_open && !self.request_window_close() {
            window_open = true;
        }
        *open = window_open;

        self.show_aot_confirmation(ctx);
        self.show_document_load_confirmation(ctx);
        self.show_window_close_confirmation(ctx, open);
        self.show_document_deletion_confirmation(ctx);
        self.show_problem_switch_confirmation(ctx);

        // Help window
        if self.show_help_window {
            egui::Window::new("Help")
                .open(&mut self.show_help_window)
                .resizable(true)
                .default_size([600.0, 400.0])
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| match self.selected_problem {
                        ProblemsEnum::BVPSimple => {
                            ui.style_mut().text_styles.insert(
                                egui::TextStyle::Body,
                                egui::FontId::new(14.0, egui::FontFamily::Proportional),
                            );
                            ui.colored_label(egui::Color32::BLACK, REACTOR_ENG_HELPER);
                        }
                        ProblemsEnum::None => {
                            ui.style_mut().text_styles.insert(
                                egui::TextStyle::Body,
                                egui::FontId::new(14.0, egui::FontFamily::Proportional),
                            );
                            ui.colored_label(
                                egui::Color32::BLACK,
                                "No help available for this problem type.",
                            );
                        }
                    });
                });
        }

        // Show plot window if available
        if let Some(plot_window) = &mut self.plot_window {
            plot_window.show(ctx);
            // Clean up when plot window is closed
            if !plot_window.visible {
                self.plot_window = None;
            }
        }
    }
}

/// Comprehensive help information for all sections and fields in reactor BVP configuration.
///
/// This map provides detailed explanations for every section and field in the SIMPLE_BVP_TEMPLATE,
/// helping users understand what each parameter means, what values it accepts, and how it affects
/// the simulation. The key format is "section_name.field_name" or just "section_name" for
/// section-level descriptions.
pub const FIELD_HELP_MAP_BVP: &[(&str, &str)] = &[
    // Section descriptions
    ("initial_guess", "Initial guess values for the numerical solver. Controls starting values for iterative solution process."),
    ("process_conditions", "Main physical and chemical parameters defining the reactor problem setup."),
    ("boundary_condition", "Initial values for all variables at the reactor inlet."),
    ("diffusion_coefficients", "Molecular diffusion coefficients for each chemical substance [m²/s]."),
    ("reactions", "Chemical reactions with Arrhenius kinetic parameters [A, n, E, Q]."),
    ("solver_settings", "Structured solver backend controls. Residuals/Jacobian assembly is exclusive: choose Lambdify or AOT, then refine matrix and symbolic backends with explicit dropdowns."),
    ("bounds", "Variable bounds for solver convergence. Format: lower_bound, upper_bound for each variable type."),
    ("rel_tolerance", "Relative convergence tolerances for each variable type in the numerical solver."),
    ("strategy_params", "Advanced solver strategy parameters for damped Newton-Raphson method."),
    ("adaptive_strategy", "Adaptive grid refinement settings for automatic mesh improvement."),
    ("grid_refinement", "Grid refinement method and parameters for mesh adaptation algorithms."),
    ("postprocessing", "Output and visualization options for results processing and data export."),
    // initial_guess fields
    ("initial_guess.universal", "Universal initial guess value applied to all variables and mesh points. Type: Float. Typical: 1e-2"),
    ("initial_guess.C", "Initial guess for concentration variables. Type: Float. Range: 0.0-1.0"),
    ("initial_guess.J", "Initial guess for flux variables. Type: Float."),
    ("initial_guess.Teta", "Initial guess for dimensionless temperature variables. Type: Float."),
    ("initial_guess.q", "Initial guess for heat flux variables. Type: Float."),
    // process_conditions fields
    ("process_conditions.problem_name", "Optional problem identifier for documentation. Type: Optional<String>. Example: Some(\"HMX_decomposition\")"),
    ("process_conditions.problem_description", "Optional detailed problem description. Type: Optional<String>. Example: Some(\"HMX thermal decomposition study\")"),
    ("process_conditions.substances", "List of chemical substances in the system. Type: Vec<String>. Example: H2O, CO2, N2"),
    ("process_conditions.t0", "Start position in dimensionless coordinates. Type: Float. Recommended: 0.0"),
    ("process_conditions.t_end", "End position in dimensionless coordinates. Type: Float. Recommended: 1.0"),
    ("process_conditions.n_steps", "Initial number of grid points for spatial discretization. Type: Usize. Range: 10-1000. Higher = more accurate but slower.
    \n if enabled grid refinement mesh points will be added so recommended number is 15-50"),
    ("process_conditions.arg", "Independent variable name (spatial coordinate). Type: String. Typical: x, z, t"),
    ("process_conditions.Tm", "Maximum temperature in reaction zone [K]. All parameters like ro, D, etc. are calculated at this
    temperature. Type: Float. Range: 300-3000K."),
    ("process_conditions.L", "Characteristic length scale [m]. Type: Float. Range: 1e-6 to 1e-4. Should span complete reaction zone"),
    ("process_conditions.dT", "Temperature difference for scaling [K]. Type: Float. Range: 100-2000K. Usually equals or slightly less then 
    \n intitial temperature (temperature from boundary condition)"),
    ("process_conditions.T_scale", "Temperature scaling factor [K]. Type: Float. Usually equals to Tm - dT. Used in dimensionless equations"),
    ("process_conditions.P", "System pressure [Pa]. Type: Float. Range: 1e3-1e8 Pa. Affects reaction rates and transport"),
    ("process_conditions.Cp", "Specific heat capacity [J/kg/K]. Type: Float. Material property for energy balance"),
    ("process_conditions.Lambda", "Thermal conductivity [W/m/K]. Type: Float. Controls heat conduction"),
    ("process_conditions.m", "Mass velocity [kg/m²/s]. Type: Float. Range: 1e-6 to 1e3. Flow rate per unit area"),
    ("process_conditions.M", "Average molar mass [kg/mol]. Type: Float.  Used in gas law calculations"),
    ("process_conditions.thermal_effects", "Heat of reaction for each reaction [J/mol]. Type: Vec<Float>. Positive=exothermic, Negative=endothermic"),
    ("process_conditions.groups", "Enable custom atomic group definitions. Type: Boolean. true=use custom groups, false=standard elements"),
    // boundary_condition fields
    ("boundary_condition.T", "Initial/boundary temperature [K]. Type: Float. Range: 200-2000K. Starting temperature of the system"),
    // solver_settings fields
    ("solver_settings.scheme", "Spatial discretization scheme. Type: String. Options: forward, trapezoid. forward remains the lightest-weight route."),
    ("solver_settings.method", "Linear algebra backend. Type: String. Options: Banded, Sparse. This is the canonical matrix-storage selector."),
    ("solver_settings.strategy", "Convergence strategy. Type: String. Options: Damped, Frozen, Naive. Damped remains the robust default."),
    ("solver_settings.linear_sys_method", "Optional linear-system override. Type: Optional<String>. Leave None for the default matrix-backend policy."),
    ("solver_settings.abs_tolerance", "Absolute convergence tolerance. Type: Float. Range: 1e-12 to 1e-3. Smaller=more accurate"),
    ("solver_settings.max_iterations", "Maximum solver iterations. Type: Usize. Range: 10-1000. Higher=more attempts to converge"),
    ("solver_settings.loglevel", "Logging verbosity level. Type: Optional<String>. Options: Some(info), Some(debug), None"),
    ("solver_settings.dont_save_log", "Disable log file saving. Type: Boolean. true=no log files, false=save logs"),
    ("solver_settings.generated_backend", "High-level solver preset used by the RustedSciThe backend. Type: String. Default: banded_lambdify."),
    ("solver_settings.matrix_backend", "Compatibility alias for solver_settings.method. The structured GUI synchronizes this field automatically."),
    ("solver_settings.backend_policy", "Derived execution policy from the visible backend selector. KiThe writes lambdify_only or aot_only automatically."),
    ("solver_settings.symbolic_backend", "Symbolic assembly backend. Type: String. Options: AtomView, ExprLegacy. AtomView is the default."),
    ("solver_settings.aot_codegen_backend", "AOT code generation backend. Type: String. Options: C, Rust, Zig. Only visible when the execution backend is AOT."),
    ("solver_settings.aot_c_compiler", "C compiler choice, visible only for the C codegen backend. Type: String. Options: tcc, gcc."),
    ("solver_settings.aot_build_policy", "AOT build lifecycle policy. Type: String. Options: use_if_available, build_if_missing, require_prebuilt, rebuild_always."),
    ("solver_settings.aot_build_profile", "AOT build profile. Type: String. Options: release, debug."),
    ("solver_settings.aot_compile_preset", "AOT compile preset. Type: String. Options: production, fast_build, dev_fastest."),
    ("solver_settings.aot_execution_policy", "Runtime execution policy for compiled callbacks. Type: String. Options: auto, sequential, parallel."),
    ("solver_settings.banded_linear_solver", "Banded linear solver policy. Type: String. Options: auto, faithful, block_tridiagonal, block_tridiagonal_consistent, faer_sparse, general_partial_pivot."),
    ("solver_settings.refinement_steps", "Iterative refinement steps for the banded linear solver. Type: Usize. Default: 5. Set 0 to disable extra refinement."),
    // bounds fields
    ("bounds.C", "Concentration variable bounds [lower, upper]. Type: Vec<Float>. Example: -10.0, 10.0. Prevents unphysical values"),
    ("bounds.J", "Flux variable bounds [lower, upper]. Type: Vec<Float>. Example: -1e20, 1e20. Wide bounds for numerical stability"),
    ("bounds.Teta", "Temperature variable bounds [lower, upper]. Type: Vec<Float>. Example: -100.0, 100.0. Dimensionless temperature limits"),
    ("bounds.q", "Heat flux variable bounds [lower, upper]. Type: Vec<Float>. Example: -1e20, 1e20. Energy flux constraints"),
    // rel_tolerance fields
    ("rel_tolerance.C", "Relative tolerance for concentration variables. Type: Float. Range: 1e-12 to 1e-3. Controls convergence precision"),
    ("rel_tolerance.J", "Relative tolerance for flux variables. Type: Float. Range: 1e-12 to 1e-3. Smaller=stricter convergence"),
    ("rel_tolerance.Teta", "Relative tolerance for temperature variables. Type: Float. Range: 1e-12 to 1e-3. Temperature convergence criterion"),
    ("rel_tolerance.q", "Relative tolerance for heat flux variables. Type: Float. Range: 1e-12 to 1e-3. Energy balance precision"),
    // strategy_params fields
    ("strategy_params.max_jac", "Maximum Jacobian matrix updates per iteration. Type: Optional<Integer>. Range: 1-10. Higher=more accurate but slower"),
    ("strategy_params.max_damp_iter", "Maximum damping iterations for convergence. Type: Optional<Integer>. Range: 1-50. Controls step size reduction"),
    ("strategy_params.damp_factor", "Damping factor for step size control. Type: Optional<Float>. Range: 0.1-0.9. Smaller=more conservative steps"),
    // adaptive_strategy fields
    ("adaptive_strategy.version", "Grid refinement contract version. KiThe currently supports and seeds version 1 automatically."),
    ("adaptive_strategy.max_refinements", "Maximum number of refinement iterations. Type: Integer. Range: 1-10. More=finer final mesh"),
    // grid_refinement fields
    ("grid_refinement.grcar_smooke", "Grcar-Smooke refinement parameters [tol1, tol2, factor]. Type: Vec<Float>. Controls mesh adaptation sensitivity"),
    ("grid_refinement.pearson", "Pearson refinement parameters [tol, factor]. Type: Vec<Float>. Alternative mesh refinement method"),
    ("grid_refinement.twopnt", "Two-point refinement parameters [tol1, tol2, factor]. Type: Vec<Float>. Boundary layer focused refinement"),
    ("grid_refinement.easy", "Simple refinement parameter [factor]. Type: Vec<Float>. Basic uniform refinement"),
    ("grid_refinement.double_points", "Double points refinement (no parameters). Type: Vec<Float>. Simple point doubling method"),
    // postprocessing fields
    ("postprocessing.gnuplot", "Generate gnuplot visualization files. Type: Boolean. true=create plots, false=no plotting"),
    ("postprocessing.save_to_csv", "Export results to CSV format. Type: Boolean. true=save CSV files, false=no CSV export"),
    ("postprocessing.filename", "Output filename prefix for saved files. Type: String. Used for all output file naming"),
    ("postprocessing.plot", "Generate native Rust plots. Type: Boolean. true=create internal plots, false=no native plotting"),
    ("postprocessing.save", "Save results to text files. Type: Boolean. true=save text output, false=no text files"),
    ("postprocessing.return_to_dimension", "Convert results back to dimensional units. Type: Boolean. true=dimensional output, false=dimensionless"),
    ("postprocessing.no_plots_in_terminal", "Disable terminal plot output. Type: Boolean. true=no terminal plots, false=show in terminal"),
];
