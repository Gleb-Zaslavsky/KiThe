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
use crate::ReactorsBVP::task_parser_reactor_BVP::SIMPLE_BVP_TEMPLATE;
use crate::cli::reactor_help::REACTOR_ENG_HELPER;
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::Utils::task_parser::{DocumentMap, DocumentParser, Value};
use eframe::egui;
use std::collections::HashMap;
/*
this is value enum from RustedSciThe::Utils::task_parser
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
/// * `Usize` - Drag value widget with automatic type conversion
/// * `Vector` - Multiple drag widgets for each vector element
/// * `Boolean` - Checkbox widget
/// * `Optional` - Collapsible section for nested values
fn render_value(ui: &mut egui::Ui, value: &mut Value, field_name: &str) {
    match value {
        Value::String(s) => {
            if ui.text_edit_singleline(s).changed() {
                println!("Field '{}' changed to: {}", field_name, s);
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
                println!("Field '{}' changed to: {}", field_name, f);
            }
        }
        Value::Integer(i) => {
            if ui.add(egui::DragValue::new(i).speed(1)).changed() {
                println!("Field '{}' changed to: {}", field_name, i);
            }
        }
        Value::Usize(u) => {
            let mut tmp = *u as i64;
            if ui.add(egui::DragValue::new(&mut tmp).speed(1)).changed() {
                *u = tmp as usize;
                println!("Field '{}' changed to: {}", field_name, u);
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
                    println!("Field '{}[{}]' changed to: {}", field_name, i, v);
                }
            }
        }
        Value::Boolean(b) => {
            if ui.checkbox(b, "").changed() {
                println!("Field '{}' changed to: {}", field_name, b);
            }
        }
        Value::Optional(opt) => {
            ui.horizontal(|ui| {
                let mut is_some = opt.is_some();
                if ui.checkbox(&mut is_some, "").changed() {
                    if is_some && opt.is_none() {
                        *opt = Some(Box::new(Value::Float(0.0)));
                        println!("Field '{}' changed to Some", field_name);
                    } else if !is_some && opt.is_some() {
                        *opt = None;
                        println!("Field '{}' changed to None", field_name);
                    }
                }

                if let Some(inner) = opt {
                    ui.label("Some:");
                    render_value(ui, inner, field_name);
                } else {
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
) -> (bool, bool) {
    let mut add_field = false;
    let mut delete_section = false;

    ui.horizontal(|ui| {
        let section_header = egui::CollapsingHeader::new(name)
            .default_open(false)
            .show(ui, |ui| {
                let mut fields_to_delete = Vec::new();

                for (field_name, maybe_values) in section.iter_mut() {
                    ui.horizontal(|ui| {
                        let field_label = ui.label(field_name);

                        // Update help info on hover
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
                            // Special handling for Vec<String> fields like substances
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

                                if ui.text_edit_singleline(&mut combined_text).changed() {
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
                                    println!(
                                        "Field '{}' changed to: {}",
                                        field_name, combined_text
                                    );
                                }
                            } else {
                                for v in values.iter_mut() {
                                    render_value(ui, v, field_name);
                                }
                            }
                        } else {
                            ui.label("None");
                        }

                        if ui.button("❌").clicked() {
                            fields_to_delete.push(field_name.clone());
                        }
                    });
                }

                // Delete fields outside the iterator
                for field_name in fields_to_delete {
                    section.remove(&field_name);
                    println!("Deleted field: {}", field_name);
                }

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

    (add_field, delete_section)
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
    pub show_help_window: bool,
    pub current_help_info: String,
    pub plot_window: Option<PlotWindow>,
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
        for (section_name, section_data) in &self.document {
            result.push_str(section_name);
            result.push('\n');
            for (field_name, maybe_values) in section_data {
                if let Some(values) = maybe_values {
                    for value in values {
                        match value {
                            Value::String(s) => {
                                result.push_str(&format!("{}: {}\n", field_name, s))
                            }
                            Value::Float(f) => result.push_str(&format!("{}: {}\n", field_name, f)),
                            Value::Integer(i) => {
                                result.push_str(&format!("{}: {}\n", field_name, i))
                            }
                            Value::Usize(u) => result.push_str(&format!("{}: {}\n", field_name, u)),
                            Value::Boolean(b) => {
                                result.push_str(&format!("{}: {}\n", field_name, b))
                            }
                            Value::Vector(vec) => {
                                let vec_str = vec
                                    .iter()
                                    .map(|v| v.to_string())
                                    .collect::<Vec<_>>()
                                    .join(", ");
                                result.push_str(&format!("{}: [{}]\n", field_name, vec_str));
                            }
                            Value::Optional(opt) => {
                                if let Some(inner) = opt {
                                    match inner.as_ref() {
                                        Value::String(s) => result
                                            .push_str(&format!("{}: Some({})\n", field_name, s)),
                                        _ => {
                                            result.push_str(&format!("{}: Some(...)\n", field_name))
                                        }
                                    }
                                } else {
                                    result.push_str(&format!("{}: None\n", field_name));
                                }
                            }
                        }
                    }
                } else {
                    result.push_str(&format!("{}: None\n", field_name));
                }
            }
        }
        result
    }

    /// Saves the current document to the specified file path.
    ///
    /// This method converts the document to string format and writes it to disk.
    /// It provides console feedback about the success or failure of the operation.
    ///
    /// # Arguments
    ///
    /// * `path` - The file path where the document should be saved
    ///
    /// # Error Handling
    ///
    /// File I/O errors are caught and logged to the console rather than
    /// propagated, maintaining GUI stability.
    pub fn save_document(&self, path: std::path::PathBuf) {
        let content = self.document_to_string();
        match std::fs::write(&path, content) {
            Ok(_) => println!("Successfully saved to: {:?}", path),
            Err(e) => println!("Error saving file: {}", e),
        }
    }

    /// Executes the appropriate calculation based on the selected problem type.
    ///
    /// This method dispatches to the correct solver based on the current problem
    /// type selection. It passes the current document state to the solver,
    /// allowing users to run calculations with their configured parameters.
    ///
    /// # Problem Type Dispatch
    ///
    /// * `BVPSimple` - Creates a SimpleReactorTask and calls solve_from_parsed()
    /// * `None` - Shows a message that no calculation is available
    ///
    /// # Error Handling
    ///
    /// Solver errors are handled internally by the respective solver implementations.
    /// This method focuses on dispatch logic rather than error propagation.
    fn run_calculation(&mut self) {
        match self.selected_problem {
            ProblemsEnum::BVPSimple => {
                println!("Starting BVP Simple calculation...");
                let mut reactor = SimpleReactorTask::new();
                let _ = reactor.solve_from_parsed(self.document.clone());
                let postproc = &self.document.get("postprocessing").cloned().unwrap();
                let gui_plot = if let Some(gui_plot_) = postproc.get("gui_plot") {
                    gui_plot_.clone().unwrap()[0]
                        .as_boolean()
                        .expect("Failed to get gui_plot as bool")
                } else {
                    true
                };

                if gui_plot {
                    let y = reactor.solver.solution.clone().unwrap();
                    let x_mesh = reactor.solver.x_mesh.clone().unwrap();
                    let arg = "x".to_owned();
                    let values = reactor.solver.unknowns.clone();

                    let t_result = x_mesh;
                    let y_result =
                        nalgebra::DMatrix::from_column_slice(y.nrows(), y.ncols(), y.as_slice());

                    self.plot_window = Some(PlotWindow::new(arg, values, t_result, y_result));
                }

                return;
            }
            ProblemsEnum::None => {
                println!("No calculation available for this problem type.");
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
                    Ok(_) => println!("Parsed successfully"),
                    Err(e) => println!("Error parsing template: {}", e),
                }
                parser.get_result().cloned().unwrap_or_else(HashMap::new)
            }
            None => HashMap::new(),
        };

        Self {
            document,
            new_section_name: String::new(),
            new_field_names: HashMap::new(),
            new_field_values: HashMap::new(),
            selected_problem: problem,
            current_file_path: None,
            show_help_window: false,
            current_help_info: String::new(),
            plot_window: None,
        }
    }

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
        egui::Window::new("Gas-phase Combustion/Steady State Plug Flow")
            .open(open)
            .default_size(ctx.screen_rect().size())
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
                                    self.save_document(path.clone());
                                } else {
                                    if let Some(path) = rfd::FileDialog::new()
                                        .add_filter("text", &["txt"])
                                        .save_file()
                                    {
                                        self.current_file_path = Some(path.clone());
                                        self.save_document(path);
                                    }
                                }
                            }
                            if ui.button("💾 Save As...").clicked() {
                                if let Some(path) = rfd::FileDialog::new()
                                    .add_filter("text", &["txt"])
                                    .save_file()
                                {
                                    self.current_file_path = Some(path.clone());
                                    self.save_document(path);
                                }
                            }
                            if ui.button("❓ Help").clicked() {
                                self.show_help_window = true;
                            }
                        });

                    ui.separator();

                    // Problem selection dropdown
                    ui.label("Problem Type:");
                    let old_problem = self.selected_problem.clone();
                    egui::ComboBox::from_label("")
                        .selected_text(format!("{}", self.selected_problem))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.selected_problem,
                                ProblemsEnum::None,
                                "None",
                            );
                            ui.selectable_value(
                                &mut self.selected_problem,
                                ProblemsEnum::BVPSimple,
                                "BVP Simple",
                            );
                        });

                    // If problem changed, reload document
                    if old_problem as u8 != self.selected_problem.clone() as u8 {
                        let new_app = Self::new_with_problem(self.selected_problem.clone());
                        self.document = new_app.document;
                        println!("Switched to problem type: {}", self.selected_problem);
                    }

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("📁 Read Task").clicked() {
                            if let Some(path) = rfd::FileDialog::new()
                                .add_filter("text", &["txt"])
                                .add_filter("all", &["*"])
                                .pick_file()
                            {
                                println!("Selected file: {:?}", path);
                                match std::fs::read_to_string(&path) {
                                    Ok(content) => {
                                        let mut parser = DocumentParser::new(content);
                                        match parser.parse_document() {
                                            Ok(_) => {
                                                if let Some(parsed_doc) = parser.get_result() {
                                                    self.document = parsed_doc.clone();
                                                    self.current_file_path = Some(path.clone());
                                                    println!(
                                                        "Successfully loaded and parsed file: {:?}",
                                                        path
                                                    );
                                                } else {
                                                    println!("Error: Parser returned no result");
                                                }
                                            }
                                            Err(e) => println!("Error parsing file: {}", e),
                                        }
                                    }
                                    Err(e) => println!("Error reading file: {}", e),
                                }
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

                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("ℹ️");
                    ui.style_mut().text_styles.insert(
                        egui::TextStyle::Body,
                        egui::FontId::new(16.0, egui::FontFamily::Proportional),
                    );
                    let help_width = ui.available_width() - 30.0;
                    let display_text = if self.current_help_info.is_empty() {
                        "Hover over fields and sections for help information"
                    } else {
                        &self.current_help_info
                    };
                    ui.add_sized(
                        [help_width, 80.0],
                        egui::Label::new(
                            egui::RichText::new(display_text)
                                .color(egui::Color32::BLACK)
                                .size(16.0),
                        )
                        .wrap(),
                    );
                });
                ui.add_space(20.0);

                // Render document sections
                egui::ScrollArea::vertical().show(ui, |ui| {
                    let mut sections_to_process = Vec::new();
                    let mut sections_to_delete = Vec::new();

                    for (section_name, section_data) in self.document.iter_mut() {
                        let field_name = self
                            .new_field_names
                            .entry(section_name.clone())
                            .or_insert_with(String::new);
                        let field_value = self
                            .new_field_values
                            .entry(section_name.clone())
                            .or_insert_with(String::new);

                        let (add_field, delete_section) = render_section(
                            ui,
                            section_name,
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
                    }

                    // Delete sections outside the iterator
                    for section_name in sections_to_delete {
                        self.document.remove(&section_name);
                        println!("Deleted section: {}", section_name);
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
                                section.insert(field_name.clone(), Some(vec![value]));
                                println!(
                                    "Added field '{}' to section '{}' with value: {}",
                                    field_name, section_name, field_value
                                );
                            }

                            self.new_field_names
                                .insert(section_name.clone(), String::new());
                            self.new_field_values.insert(section_name, String::new());
                        }
                    }

                    ui.separator();
                    ui.horizontal(|ui| {
                        ui.label("New section name:");
                        ui.text_edit_singleline(&mut self.new_section_name);
                        if ui.button("Create Section").clicked()
                            && !self.new_section_name.is_empty()
                        {
                            self.document
                                .insert(self.new_section_name.clone(), HashMap::new());
                            println!("Created new section: {}", self.new_section_name);
                            self.new_section_name.clear();
                        }
                    });
                });

                // Big RUN CALCULATION button at the bottom
                ui.separator();
                ui.add_space(10.0);
                ui.with_layout(egui::Layout::top_down(egui::Align::Center), |ui| {
                    let button = egui::Button::new("🚀 RUN CALCULATION!")
                        .min_size(egui::Vec2::new(200.0, 50.0))
                        .fill(egui::Color32::from_rgb(160, 160, 160));

                    if ui.add_sized([200.0, 50.0], button).clicked() {
                        self.run_calculation();
                    }
                });
            });

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
    ("solver_settings", "Numerical solver configuration including discretization scheme and solution method."),
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
    ("solver_settings.scheme", "Spatial discretization scheme. Type: String. Options: forward, backward, central. forward=upwind stable"),
    ("solver_settings.method", "Linear system solution method. Type: String. Options: Sparse, Dense. Sparse=efficient for large systems"),
    ("solver_settings.strategy", "Convergence strategy. Type: String. Options: Damped, Newton. Damped=more robust convergence"),
    ("solver_settings.linear_sys_method", "Linear solver algorithm. Type: Optional<String>. None=automatic selection"),
    ("solver_settings.abs_tolerance", "Absolute convergence tolerance. Type: Float. Range: 1e-12 to 1e-3. Smaller=more accurate"),
    ("solver_settings.max_iterations", "Maximum solver iterations. Type: Integer. Range: 10-1000. Higher=more attempts to converge"),
    ("solver_settings.loglevel", "Logging verbosity level. Type: Optional<String>. Options: Some(info), Some(debug), None"),
    ("solver_settings.dont_save_logs", "Disable log file saving. Type: Boolean. true=no log files, false=save logs"),
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
    ("adaptive_strategy.version", "Grid refinement algorithm version. Type: Integer. Options: 1, 2. Different refinement strategies"),
    ("adaptive_strategy.max_refinements", "Maximum number of refinement iterations. Type: Integer. Range: 1-10. More=finer final mesh"),
    // grid_refinement fields
    ("grid_refinement.grcarsmooke", "Grcar-Smooke refinement parameters [tol1, tol2, factor]. Type: Vec<Float>. Controls mesh adaptation sensitivity"),
    ("grid_refinement.pearson", "Pearson refinement parameters [tol, factor]. Type: Vec<Float>. Alternative mesh refinement method"),
    ("grid_refinement.twopnt", "Two-point refinement parameters [tol1, tol2, factor]. Type: Vec<Float>. Boundary layer focused refinement"),
    ("grid_refinement.easy", "Simple refinement parameter [factor]. Type: Vec<Float>. Basic uniform refinement"),
    ("grid_refinement.doubleoints", "Double points refinement (no parameters). Type: Vec<Float>. Simple point doubling method"),
    // postprocessing fields
    ("postprocessing.gnuplot", "Generate gnuplot visualization files. Type: Boolean. true=create plots, false=no plotting"),
    ("postprocessing.save_to_csv", "Export results to CSV format. Type: Boolean. true=save CSV files, false=no CSV export"),
    ("postprocessing.filename", "Output filename prefix for saved files. Type: String. Used for all output file naming"),
    ("postprocessing.plot", "Generate native Rust plots. Type: Boolean. true=create internal plots, false=no native plotting"),
    ("postprocessing.save", "Save results to text files. Type: Boolean. true=save text output, false=no text files"),
    ("postprocessing.return_to_dimension", "Convert results back to dimensional units. Type: Boolean. true=dimensional output, false=dimensionless"),
    ("postprocessing.no_plots_in_terminal", "Disable terminal plot output. Type: Boolean. true=no terminal plots, false=show in terminal"),
];
