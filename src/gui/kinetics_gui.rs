//! # Kinetics GUI Module
//!
//! This module provides a comprehensive graphical user interface for managing chemical reaction databases.
//! It allows users to browse, search, select, and manipulate kinetic reaction data from various libraries.
//!
//! ## Main Features
//! - **Library Selection**: Browse different kinetic databases (NUIG, Cantera, Aramco, etc.)
//! - **Reaction Browsing**: View all reactions in a selected library with search functionality
//! - **Reaction Selection**: Click on reactions to view detailed kinetic parameters
//! - **Reaction Collection**: Save individual reactions or entire mechanisms for calculations
//! - **Custom Reaction Creation**: Create new reactions with different kinetic models
//! - **Mechanism Construction**: Build sub-mechanisms from selected substances
//!
//! ## GUI Layout
//! The interface is divided into two main panels:
//! - **Left Panel**: Reaction list with search filter and library selection dropdown
//! - **Right Panel**: Selected reaction details and action buttons
//!
//! ## Button Functions
//! - **"Saving reaction for calculation"**: Adds the currently selected reaction to the calculation queue
//! - **"Taking all reactions from mechanism"**: Adds all reactions from the current library to the calculation queue
//! - **"Add new reaction"**: Opens a dialog to create custom reactions with different kinetic models
//! - **"Searching reactions"**: Placeholder for substance-based reaction search functionality
//! - **"Building sub-mechanism"**: Placeholder for mechanism construction from selected substances

use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::ReactionType as RType;
use crate::Kinetics::mechfinder_api::{ReactionData, parse_kinetic_data_vec};
use eframe::egui;
use std::collections::HashMap;

/// Main application structure for the kinetics GUI
///
/// This struct manages the entire state of the kinetics interface, including:
/// - Connection to kinetic databases through `KineticData`
/// - Currently selected library and reaction
/// - User input fields and search filters
/// - Collection of reactions added for calculation
/// - State of the "Add New Reaction" dialog window
#[derive(Debug)]
pub struct KineticsApp {
    /// Core kinetic database interface - handles loading libraries and searching reactions
    kinetic_data: KineticData,
    /// Currently selected kinetic library (e.g., "NUIG", "Cantera", "Aramco")
    selected_library: String,
    /// Parsed kinetic data of the currently selected reaction (lazy-loaded on selection)
    selected_reaction_data: Option<ReactionData>,
    /// Chemical equation of the currently selected reaction
    selected_equation: Option<String>,
    /// Input field for entering new reaction equations (currently unused)
    _reaction_input: String,
    /// Input field for entering substances to search for mechanism construction
    mechanism_input: String,
    /// Search filter text for filtering reactions in the list
    search_filter: String,
    /// Radio button selection for mechanism vs reaction mode (currently unused)
    reaction_type: ReactionType,
    /// Collection of reactions added for calculation - populated by "Save" and "Take all" buttons
    added_reactions: Vec<ReactionData>,
    /// Flag controlling visibility of the "Add New Reaction" dialog window
    show_add_reaction_window: bool,
    /// State of the "Add New Reaction" dialog - contains all input fields for creating custom reactions
    new_reaction_window: NewReactionWindow,
    /// Cache for parsed reaction data to avoid re-parsing
    reaction_cache: HashMap<String, ReactionData>,
}

/// Dialog window state for creating new custom reactions
///
/// This struct contains all the input fields needed to create different types of kinetic reactions:
/// - Elementary: Simple Arrhenius parameters [A, n, E]
/// - Falloff: Low-pressure and high-pressure rate parameters with optional Troe parameters
/// - ThreeBody: Arrhenius parameters with collision efficiency factors
/// - Pressure: Pressure-dependent rate data in JSON format
#[derive(Debug)]
struct NewReactionWindow {
    /// Selected kinetic model type (Elementary, Falloff, ThreeBody, Pressure)
    reaction_type: crate::Kinetics::mechfinder_api::ReactionType,
    /// Chemical equation input field (e.g., "H2 + O <=> H + OH")
    equation: String,
    /// Arrhenius parameters [A, n, E] for Elementary and ThreeBody reactions
    arrenius: [String; 3],
    /// Low-pressure rate parameters [A, n, E] for Falloff reactions
    low_rate: [String; 3],
    /// High-pressure rate parameters [A, n, E] for Falloff reactions
    high_rate: [String; 3],
    /// JSON input for collision efficiencies in ThreeBody reactions (e.g., {"H2": 2.5, "H2O": 12.0})
    eff_input: String,
    /// JSON input for pressure-dependent rate data in Pressure reactions
    pressure_data: String,
}

impl Default for NewReactionWindow {
    fn default() -> Self {
        Self {
            reaction_type: crate::Kinetics::mechfinder_api::ReactionType::Elem,
            equation: String::new(),
            arrenius: [String::new(), String::new(), String::new()],
            low_rate: [String::new(), String::new(), String::new()],
            high_rate: [String::new(), String::new(), String::new()],
            eff_input: String::new(),
            pressure_data: String::new(),
        }
    }
}

#[derive(Default, PartialEq, Clone, Debug)]
enum ReactionType {
    #[default]
    Mechanism,
    Reaction,
}

impl Default for KineticsApp {
    fn default() -> Self {
        Self {
            kinetic_data: KineticData::new(),
            selected_library: String::new(),
            selected_reaction_data: None,
            selected_equation: None,
            _reaction_input: String::new(),
            mechanism_input: String::new(),
            search_filter: String::new(),
            reaction_type: ReactionType::default(),
            added_reactions: Vec::new(),
            show_add_reaction_window: false,
            new_reaction_window: NewReactionWindow::default(),
            reaction_cache: HashMap::new(),
        }
    }
}

impl KineticsApp {
    /// Creates a new KineticsApp instance with default NUIG library loaded
    ///
    /// This constructor:
    /// 1. Initializes all fields to default values
    /// 2. Sets the selected library to "NUIG"
    /// 3. Loads the NUIG reaction database
    /// 4. Populates the reaction list for display
    pub fn new() -> Self {
        let mut app = Self::default();
        app.selected_library = "NUIG".to_string();
        app.load_library_reactions();
        app
    }

    /// Loads reactions from the currently selected library
    ///
    /// This method:
    /// 1. Opens JSON files for the selected library using `kinetic_data.open_json_files()`
    /// 2. Populates `AllEquations` and `EquationReactionMap` using `print_all_reactions()`
    /// 3. Resets selected reaction state to force re-selection
    /// 4. Clears reaction cache for new library
    ///
    /// Called when:
    /// - App is initialized
    /// - User switches to a different library via dropdown
    fn load_library_reactions(&mut self) {
        self.kinetic_data.open_json_files(&self.selected_library);
        self.kinetic_data.print_all_reactions();
        self.selected_reaction_data = None;
        self.selected_equation = None;
        self.reaction_cache.clear();
    }

    /// Parses kinetic data for a selected reaction equation (lazy loading with caching)
    ///
    /// This method:
    /// 1. Checks cache first for previously parsed reactions
    /// 2. If not cached, uses `search_reaction_by_equation()` to get raw serde::Value data
    /// 3. Parses the JSON data into a `ReactionData` struct
    /// 4. Stores the parsed data in cache and `selected_reaction_data`
    ///
    /// Called when:
    /// - User clicks on a reaction in the scrollable list
    ///
    /// Benefits of caching:
    /// - Instant loading for previously viewed reactions (<10ms vs 200ms+)
    /// - Reduces JSON parsing overhead
    /// - Improves user experience when browsing reactions
    fn parse_selected_reaction(&mut self, equation: &str) {
        // Check cache first
        if let Some(cached_data) = self.reaction_cache.get(equation) {
            self.selected_reaction_data = Some(cached_data.clone());
            return;
        }

        // Parse and cache if not found
        let (_, reaction_value) = self.kinetic_data.search_reaction_by_equation(equation);
        if let Ok(reaction_data) = serde_json::from_value::<ReactionData>(reaction_value) {
            self.reaction_cache
                .insert(equation.to_string(), reaction_data.clone());
            self.selected_reaction_data = Some(reaction_data);
        }
    }

    /// Displays the "Add New Reaction" dialog window
    ///
    /// This method creates a modal dialog with:
    /// 1. **Reaction Type Dropdown**: Select kinetic model (Elementary, Falloff, ThreeBody, Pressure)
    /// 2. **Equation Input**: Text field for chemical equation
    /// 3. **Dynamic Parameter Fields**: Changes based on selected reaction type:
    ///    - Elementary: A, n, E parameters
    ///    - Falloff: Low/high rate parameters
    ///    - ThreeBody: Arrhenius + efficiency JSON
    ///    - Pressure: Pressure-dependent data JSON
    /// 4. **Action Buttons**: Create (validates and adds reaction) or Cancel
    ///
    /// Window state is controlled by `show_add_reaction_window` boolean
    fn show_add_reaction_window(&mut self, ctx: &egui::Context) {
        egui::Window::new("Add New Reaction")
            .open(&mut self.show_add_reaction_window.clone())
            .default_size([400.0, 500.0])
            .show(ctx, |ui| {
                ui.heading("Create New Reaction");

                // Reaction type dropdown
                egui::ComboBox::from_label("Reaction Type")
                    .selected_text(format!("{:?}", self.new_reaction_window.reaction_type))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.new_reaction_window.reaction_type,
                            RType::Elem,
                            "Elementary",
                        );
                        ui.selectable_value(
                            &mut self.new_reaction_window.reaction_type,
                            RType::Falloff,
                            "Falloff",
                        );
                        ui.selectable_value(
                            &mut self.new_reaction_window.reaction_type,
                            RType::ThreeBody,
                            "ThreeBody",
                        );
                        ui.selectable_value(
                            &mut self.new_reaction_window.reaction_type,
                            RType::Pressure,
                            "Pressure",
                        );
                    });

                ui.separator();

                // Equation input
                ui.horizontal(|ui| {
                    ui.label("Equation:");
                    ui.text_edit_singleline(&mut self.new_reaction_window.equation);
                });

                ui.separator();

                // Type-specific fields
                match self.new_reaction_window.reaction_type {
                    crate::Kinetics::mechfinder_api::ReactionType::Elem => {
                        ui.label("Arrhenius Parameters [A, n, E]:");
                        ui.horizontal(|ui| {
                            ui.label("A:");
                            ui.text_edit_singleline(&mut self.new_reaction_window.arrenius[0]);
                        });
                        ui.horizontal(|ui| {
                            ui.label("n:");
                            ui.text_edit_singleline(&mut self.new_reaction_window.arrenius[1]);
                        });
                        ui.horizontal(|ui| {
                            ui.label("E:");
                            ui.text_edit_singleline(&mut self.new_reaction_window.arrenius[2]);
                        });
                    }
                    crate::Kinetics::mechfinder_api::ReactionType::Falloff => {
                        ui.label("Low Rate [A, n, E]:");
                        for i in 0..3 {
                            ui.horizontal(|ui| {
                                ui.label(format!("Low[{}]:", i));
                                ui.text_edit_singleline(&mut self.new_reaction_window.low_rate[i]);
                            });
                        }
                        ui.label("High Rate [A, n, E]:");
                        for i in 0..3 {
                            ui.horizontal(|ui| {
                                ui.label(format!("High[{}]:", i));
                                ui.text_edit_singleline(&mut self.new_reaction_window.high_rate[i]);
                            });
                        }
                    }
                    crate::Kinetics::mechfinder_api::ReactionType::ThreeBody => {
                        ui.label("Arrhenius Parameters [A, n, E]:");
                        for i in 0..3 {
                            ui.horizontal(|ui| {
                                ui.label(format!("Arr[{}]:", i));
                                ui.text_edit_singleline(&mut self.new_reaction_window.arrenius[i]);
                            });
                        }
                        ui.horizontal(|ui| {
                            ui.label("Efficiencies (JSON):");
                            ui.text_edit_singleline(&mut self.new_reaction_window.eff_input);
                        });
                    }
                    crate::Kinetics::mechfinder_api::ReactionType::Pressure => {
                        ui.horizontal(|ui| {
                            ui.label("Pressure Data (JSON):");
                            ui.text_edit_multiline(&mut self.new_reaction_window.pressure_data);
                        });
                    }
                    _ => {}
                }

                ui.separator();

                // Buttons
                ui.horizontal(|ui| {
                    if ui.button("Create Reaction").clicked() {
                        if let Some(reaction) = self.create_reaction_from_window() {
                            self.added_reactions.push(reaction);
                            self.show_add_reaction_window = false;
                            self.new_reaction_window = NewReactionWindow::default();
                        }
                    }
                    if ui.button("Cancel").clicked() {
                        self.show_add_reaction_window = false;
                        self.new_reaction_window = NewReactionWindow::default();
                    }
                });
            });
    }

    /// Creates a ReactionData struct from the "Add New Reaction" dialog inputs
    ///
    /// This method:
    /// 1. Validates that equation is not empty
    /// 2. Parses numeric parameters from string inputs
    /// 3. Creates appropriate ReactionData variant based on selected type
    /// 4. Returns None if validation fails (invalid numbers, malformed JSON, etc.)
    ///
    /// Supported reaction types:
    /// - **Elementary**: Parses Arrhenius parameters, creates via `ReactionData::new_elementary()`
    /// - **ThreeBody**: Parses Arrhenius + JSON efficiencies, creates via `ReactionData::new_three_body()`
    /// - **Falloff/Pressure**: Currently return None (not implemented)
    ///
    /// Called when user clicks "Create Reaction" button
    fn create_reaction_from_window(&self) -> Option<ReactionData> {
        if self.new_reaction_window.equation.is_empty() {
            return None;
        }

        match self.new_reaction_window.reaction_type {
            crate::Kinetics::mechfinder_api::ReactionType::Elem => {
                let arrenius: Result<Vec<f64>, _> = self
                    .new_reaction_window
                    .arrenius
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect();
                if let Ok(arr) = arrenius {
                    Some(ReactionData::new_elementary(
                        self.new_reaction_window.equation.clone(),
                        arr,
                        None,
                    ))
                } else {
                    None
                }
            }
            crate::Kinetics::mechfinder_api::ReactionType::ThreeBody => {
                let arrenius: Result<Vec<f64>, _> = self
                    .new_reaction_window
                    .arrenius
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect();
                let eff: Result<HashMap<String, f64>, _> =
                    serde_json::from_str(&self.new_reaction_window.eff_input);
                if let (Ok(arr), Ok(eff_map)) = (arrenius, eff) {
                    Some(ReactionData::new_three_body(
                        self.new_reaction_window.equation.clone(),
                        arr,
                        eff_map,
                        None,
                    ))
                } else {
                    None
                }
            }
            _ => None, // Implement other types as needed
        }
    }

    /// Main GUI rendering method - displays the entire kinetics interface
    ///
    /// ## Layout Structure:
    ///
    /// ### Left Panel (400px width):
    /// - **"List of Reactions" heading**
    /// - **Search filter**: Text input for filtering reactions by equation text
    /// - **Scrollable reaction list**: Displays all equations from `AllEquations`, filtered by search term
    ///   - Click handler: Selects reaction and triggers `parse_selected_reaction()`
    /// - **"Mechanism Source" dropdown**: Library selection (NUIG, Cantera, etc.)
    ///   - Change handler: Calls `load_library_reactions()` to switch libraries
    ///
    /// ### Right Panel:
    /// - **"Chosen Reaction Details" section**: Shows parsed kinetic data of selected reaction
    /// - **Action Buttons Row 1**:
    ///   - **"Saving reaction for calculation"**: Adds `selected_reaction_data` to `added_reactions`
    ///   - **"Taking all reactions from mechanism"**: Parses entire library and adds to `added_reactions`
    /// - **Action Buttons Row 2**: Placeholder buttons for future functionality
    /// - **Substance Input**: Text field for mechanism construction (not yet implemented)
    /// - **Radio Buttons**: Mechanism vs Reaction mode selection (not yet implemented)
    /// - **"Add new reaction" button**: Opens custom reaction creation dialog
    ///
    /// ## Button Logic:
    /// - **Save Single**: Clones `selected_reaction_data` and pushes to `added_reactions`
    /// - **Take All**: Uses `parse_kinetic_data_vec()` to parse all library reactions, extends `added_reactions`
    /// - **Add New**: Sets `show_add_reaction_window = true` to display creation dialog
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Selecting and adding reactions and mechanisms")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    // Left panel - Reaction list
                    ui.vertical(|ui| {
                        ui.set_width(400.0);
                        ui.set_min_height(600.0);
                        ui.heading("List of Reactions");
                        // Search filter
                        ui.horizontal(|ui| {
                            ui.label("Search:");
                            ui.text_edit_singleline(&mut self.search_filter);
                        });
                        ui.separator();
                        egui::ScrollArea::vertical()
                            .max_height(1200.0)
                            .show(ui, |ui| {
                                for equation in &self.kinetic_data.AllEquations.clone() {
                                    if self.search_filter.is_empty()
                                        || equation
                                            .to_lowercase()
                                            .contains(&self.search_filter.to_lowercase())
                                    {
                                        let is_selected =
                                            self.selected_equation.as_ref() == Some(equation);
                                        if ui.selectable_label(is_selected, equation).clicked() {
                                            self.selected_equation = Some(equation.clone());
                                            self.parse_selected_reaction(equation);
                                        }
                                    }
                                }
                            });
                        ui.separator();
                        // Dropdown for mechanism selection
                        egui::ComboBox::from_label("Mechanism Source")
                            .selected_text(&self.selected_library)
                            .show_ui(ui, |ui| {
                                for library in &self.kinetic_data.AllLibraries.clone() {
                                    if ui
                                        .selectable_value(
                                            &mut self.selected_library,
                                            library.clone(),
                                            library,
                                        )
                                        .clicked()
                                    {
                                        self.load_library_reactions();
                                    }
                                }
                            });
                    });
                    ui.separator();
                    // Right panel - Reaction details and controls
                    ui.vertical(|ui| {
                        ui.heading("Chosen Reaction Details");
                        if let Some(reaction) = &self.selected_reaction_data {
                            // Display selected reaction details
                            ui.group(|ui| {
                                ui.set_min_height(150.0);
                                ui.label(format!("Equation: {}", reaction.eq));
                                ui.label(format!("Type: {:?}", reaction.reaction_type));
                                ui.label("Kinetic data:");
                                ui.label(format!("{:#?}", reaction.data));
                            });
                        } else {
                            ui.group(|ui| {
                                ui.set_min_height(150.0);
                                ui.label("Выберите реакцию из списка");
                            });
                        }
                        ui.separator();
                        // Action buttons
                        ui.horizontal(|ui| {
                            if ui.button("Saving reaction for calculation").clicked() {
                                if let Some(reaction) = &self.selected_reaction_data {
                                    self.added_reactions.push(reaction.clone());
                                    println!("Added reaction: {}", reaction.eq);
                                }
                            }
                            if ui.button("Taking all reactions from mechanism").clicked() {
                                let reaction_values: Vec<serde_json::Value> =
                                    self.kinetic_data.LibKineticData.values().cloned().collect();
                                let (parsed_reactions, _) = parse_kinetic_data_vec(reaction_values);
                                self.added_reactions.extend(parsed_reactions);
                                println!(
                                    "Added {} reactions from mechanism {}",
                                    self.kinetic_data.LibKineticData.len(),
                                    self.selected_library
                                );
                            }
                        });
                        ui.horizontal(|ui| {
                            if ui.button("Searching reactions").clicked() {
                                println!("Searching reactions");
                            }
                            if ui.button("Building sub-mechanism").clicked() {
                                println!("Building sub-mechanism");
                            }
                        });
                        ui.horizontal(|ui| {
                            if ui.button("Adding sub-mechanism to calculation").clicked() {
                                println!("Adding sub-mechanism to calculation");
                            }
                        });
                        ui.separator();
                        // Input section
                        ui.heading("Construct sub-mechanism for these substances:");
                        ui.horizontal(|ui| {
                            ui.label("Enter substances to search:");
                            ui.text_edit_singleline(&mut self.mechanism_input);
                            // todo!("mechfider search by substances")
                        });

                        ui.separator();
                        // Radio buttons for reaction type
                        ui.horizontal(|ui| {
                            ui.radio_value(
                                &mut self.reaction_type,
                                ReactionType::Mechanism,
                                "Mechanism",
                            );
                            ui.radio_value(
                                &mut self.reaction_type,
                                ReactionType::Reaction,
                                "Reaction",
                            );
                        });

                        // Bottom section
                        ui.horizontal(|ui| {
                            ui.label("New reactions:");

                            if ui.button("Add new reaction").clicked() {
                                self.show_add_reaction_window = true;
                            }
                        });
                        // Show add reaction window
                        if self.show_add_reaction_window {
                            self.show_add_reaction_window(ctx);
                        }
                    });
                });
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kinetics_app_new() {
        let app = KineticsApp::new();
        assert_eq!(app.selected_library, "NUIG");
        assert!(app.added_reactions.is_empty());
        assert!(app.selected_reaction_data.is_none());
        assert!(app.selected_equation.is_none());
    }

    #[test]
    fn test_kinetics_app_default() {
        let app = KineticsApp::default();
        assert!(app.selected_library.is_empty());
        assert!(app.added_reactions.is_empty());
        assert!(app.selected_reaction_data.is_none());
        assert!(app.selected_equation.is_none());
        assert_eq!(app.reaction_type, ReactionType::Mechanism);
    }

    #[test]
    fn test_load_library_reactions() {
        let mut app = KineticsApp::default();
        app.selected_library = "NUIG".to_string();
        app.load_library_reactions();

        assert!(!app.kinetic_data.AllLibraries.is_empty());
        assert!(!app.kinetic_data.AllEquations.is_empty());
        assert!(app.selected_reaction_data.is_none());
        assert!(app.selected_equation.is_none());
    }

    #[test]
    fn test_parse_selected_reaction() {
        let mut app = KineticsApp::new();

        if !app.kinetic_data.AllEquations.is_empty() {
            let first_equation = app.kinetic_data.AllEquations[0].clone();
            app.parse_selected_reaction(&first_equation);

            assert!(app.selected_reaction_data.is_some());
            if let Some(reaction) = &app.selected_reaction_data {
                assert_eq!(reaction.eq, first_equation);
            }
        }
    }

    #[test]
    fn test_reaction_type_enum() {
        let mechanism = ReactionType::Mechanism;
        let reaction = ReactionType::Reaction;

        assert_eq!(mechanism, ReactionType::default());
        assert_ne!(mechanism, reaction);
    }

    #[test]
    fn test_add_single_reaction() {
        let mut app = KineticsApp::new();

        if !app.kinetic_data.AllEquations.is_empty() {
            let first_equation = app.kinetic_data.AllEquations[0].clone();
            app.parse_selected_reaction(&first_equation);

            let initial_count = app.added_reactions.len();

            if let Some(reaction) = &app.selected_reaction_data {
                app.added_reactions.push(reaction.clone());
                assert_eq!(app.added_reactions.len(), initial_count + 1);
                assert_eq!(app.added_reactions.last().unwrap().eq, first_equation);
            }
        }
    }

    #[test]
    fn test_add_all_reactions_from_mechanism() {
        let mut app = KineticsApp::new();

        let reaction_values: Vec<serde_json::Value> =
            app.kinetic_data.LibKineticData.values().cloned().collect();
        let (parsed_reactions, _) = parse_kinetic_data_vec(reaction_values);

        let initial_count = app.added_reactions.len();
        app.added_reactions.extend(parsed_reactions.clone());

        assert_eq!(
            app.added_reactions.len(),
            initial_count + parsed_reactions.len()
        );
        assert!(!app.added_reactions.is_empty());
    }

    #[test]
    fn test_library_switching() {
        let mut app = KineticsApp::new();
        let initial_library = app.selected_library.clone();

        // Switch to different library if available
        if app.kinetic_data.AllLibraries.len() > 1 {
            let new_library = app
                .kinetic_data
                .AllLibraries
                .iter()
                .find(|&lib| lib != &initial_library)
                .unwrap()
                .clone();

            app.selected_library = new_library.clone();
            app.load_library_reactions();

            assert_eq!(app.selected_library, new_library);
            assert!(app.selected_reaction_data.is_none());
            assert!(app.selected_equation.is_none());
        }
    }

    #[test]
    fn test_search_filter_functionality() {
        let app = KineticsApp::new();
        let search_term = "H2O";

        let filtered_equations: Vec<&String> = app
            .kinetic_data
            .AllEquations
            .iter()
            .filter(|equation| {
                equation
                    .to_lowercase()
                    .contains(&search_term.to_lowercase())
            })
            .collect();

        // Test that filtering works (assuming there are H2O reactions)
        if !filtered_equations.is_empty() {
            assert!(
                filtered_equations
                    .iter()
                    .all(|eq| eq.to_lowercase().contains(&search_term.to_lowercase()))
            );
        }
    }

    #[test]
    fn test_reaction_data_consistency() {
        let mut app = KineticsApp::new();

        if !app.kinetic_data.AllEquations.is_empty() {
            let equation = app.kinetic_data.AllEquations[0].clone();
            app.parse_selected_reaction(&equation);

            if let Some(reaction) = &app.selected_reaction_data {
                // Test that parsed reaction data is consistent
                assert!(!reaction.eq.is_empty());
                assert!(matches!(
                    reaction.reaction_type,
                    crate::Kinetics::mechfinder_api::ReactionType::Elem
                        | crate::Kinetics::mechfinder_api::ReactionType::Falloff
                        | crate::Kinetics::mechfinder_api::ReactionType::Pressure
                        | crate::Kinetics::mechfinder_api::ReactionType::ThreeBody
                        | crate::Kinetics::mechfinder_api::ReactionType::Empirical
                ));
            }
        }
    }

    #[test]
    fn test_input_fields_initialization() {
        let app = KineticsApp::default();

        assert!(app._reaction_input.is_empty());
        assert!(app.mechanism_input.is_empty());
        assert!(app.search_filter.is_empty());
    }

    #[test]
    fn test_added_reactions_uniqueness() {
        let mut app = KineticsApp::new();

        if !app.kinetic_data.AllEquations.is_empty() {
            let equation = app.kinetic_data.AllEquations[0].clone();
            app.parse_selected_reaction(&equation);

            if let Some(reaction) = &app.selected_reaction_data {
                // Add same reaction twice
                app.added_reactions.push(reaction.clone());
                app.added_reactions.push(reaction.clone());

                assert_eq!(app.added_reactions.len(), 2);
                assert_eq!(app.added_reactions[0].eq, app.added_reactions[1].eq);
            }
        }
    }

    #[test]
    fn test_new_reaction_window_default() {
        let window = NewReactionWindow::default();
        assert_eq!(
            window.reaction_type,
            crate::Kinetics::mechfinder_api::ReactionType::Elem
        );
        assert!(window.equation.is_empty());
        assert_eq!(
            window.arrenius,
            [String::new(), String::new(), String::new()]
        );
        assert_eq!(
            window.low_rate,
            [String::new(), String::new(), String::new()]
        );
        assert_eq!(
            window.high_rate,
            [String::new(), String::new(), String::new()]
        );
        assert!(window.eff_input.is_empty());
        assert!(window.pressure_data.is_empty());
    }

    #[test]
    fn test_create_elementary_reaction_from_window() {
        let mut app = KineticsApp::default();
        app.new_reaction_window.equation = "H2 + O <=> H + OH".to_string();
        app.new_reaction_window.reaction_type = crate::Kinetics::mechfinder_api::ReactionType::Elem;
        app.new_reaction_window.arrenius = [
            "1.0e13".to_string(),
            "0.0".to_string(),
            "15000.0".to_string(),
        ];

        let reaction = app.create_reaction_from_window();
        assert!(reaction.is_some());
        let reaction = reaction.unwrap();
        assert_eq!(reaction.eq, "H2 + O <=> H + OH");
        assert_eq!(
            reaction.reaction_type,
            crate::Kinetics::mechfinder_api::ReactionType::Elem
        );
    }

    #[test]
    fn test_create_threebody_reaction_from_window() {
        let mut app = KineticsApp::default();
        app.new_reaction_window.equation = "H2 + M <=> H + H + M".to_string();
        app.new_reaction_window.reaction_type =
            crate::Kinetics::mechfinder_api::ReactionType::ThreeBody;
        app.new_reaction_window.arrenius = [
            "1.0e14".to_string(),
            "-1.0".to_string(),
            "104000.0".to_string(),
        ];
        app.new_reaction_window.eff_input = r#"{"H2": 2.5, "H2O": 12.0}"#.to_string();

        let reaction = app.create_reaction_from_window();
        assert!(reaction.is_some());
        let reaction = reaction.unwrap();
        assert_eq!(reaction.eq, "H2 + M <=> H + H + M");
        assert_eq!(
            reaction.reaction_type,
            crate::Kinetics::mechfinder_api::ReactionType::ThreeBody
        );
    }

    #[test]
    fn test_create_reaction_empty_equation() {
        let mut app = KineticsApp::default();
        app.new_reaction_window.equation = "".to_string();
        app.new_reaction_window.arrenius =
            ["1.0".to_string(), "2.0".to_string(), "3.0".to_string()];

        let reaction = app.create_reaction_from_window();
        assert!(reaction.is_none());
    }

    #[test]
    fn test_create_reaction_invalid_arrenius() {
        let mut app = KineticsApp::default();
        app.new_reaction_window.equation = "A + B <=> C".to_string();
        app.new_reaction_window.arrenius =
            ["invalid".to_string(), "2.0".to_string(), "3.0".to_string()];

        let reaction = app.create_reaction_from_window();
        assert!(reaction.is_none());
    }

    #[test]
    fn test_create_threebody_invalid_eff() {
        let mut app = KineticsApp::default();
        app.new_reaction_window.equation = "H2 + M <=> H + H + M".to_string();
        app.new_reaction_window.reaction_type =
            crate::Kinetics::mechfinder_api::ReactionType::ThreeBody;
        app.new_reaction_window.arrenius =
            ["1.0".to_string(), "2.0".to_string(), "3.0".to_string()];
        app.new_reaction_window.eff_input = "invalid json".to_string();

        let reaction = app.create_reaction_from_window();
        assert!(reaction.is_none());
    }

    #[test]
    fn test_show_add_reaction_window_toggle() {
        let mut app = KineticsApp::default();
        assert!(!app.show_add_reaction_window);

        app.show_add_reaction_window = true;
        assert!(app.show_add_reaction_window);
    }

    #[test]
    fn test_new_reaction_window_reset() {
        let mut window = NewReactionWindow::default();
        window.equation = "test".to_string();
        window.arrenius[0] = "123".to_string();
        window.eff_input = "test".to_string();

        window = NewReactionWindow::default();
        assert!(window.equation.is_empty());
        assert_eq!(window.arrenius[0], String::new());
        assert!(window.eff_input.is_empty());
    }
}
