// MODULE UNDER DEVELOPMENT
use eframe::egui;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Reaction {
    pub equation: String,
    pub coefficients: Vec<f64>,
    pub arrhenius_params: ArrheniusParams,
}

#[derive(Debug, Clone)]
pub struct ArrheniusParams {
    pub a: f64,
    pub n: f64,
    pub ea: f64,
}

#[derive(Default)]
pub struct KineticsApp {
    reactions: HashMap<String, Reaction>,
    reaction_keys: Vec<String>,
    selected_reaction_key: Option<String>,
    reaction_input: String,
    mechanism_input: String,
    search_filter: String,
    reaction_type: ReactionType,
}

#[derive(Default, PartialEq)]
enum ReactionType {
    #[default]
    Mechanism,
    Reaction,
}

impl Default for Reaction {
    fn default() -> Self {
        Self {
            equation: String::new(),
            coefficients: vec![0.0; 4],
            arrhenius_params: ArrheniusParams {
                a: 0.0,
                n: 0.0,
                ea: 0.0,
            },
        }
    }
}

impl KineticsApp {
    pub fn new() -> Self {
        let mut app = Self::default();
        app.load_sample_reactions();
        app
    }

    fn load_sample_reactions(&mut self) {
        let sample_reactions = vec![
            (
                "H2 + O <=> H + OH",
                ArrheniusParams {
                    a: 3.87e4,
                    n: 2.7,
                    ea: 6260.0,
                },
            ),
            (
                "HO2 + O <=> O2 + OH",
                ArrheniusParams {
                    a: 2.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "H2O2 + O <=> HO2 + OH",
                ArrheniusParams {
                    a: 9.63e6,
                    n: 2.0,
                    ea: 4000.0,
                },
            ),
            (
                "CH + O <=> CO + H",
                ArrheniusParams {
                    a: 5.7e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH2 + O <=> H + HCO",
                ArrheniusParams {
                    a: 8.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH2(S) + O <=> CO + H2",
                ArrheniusParams {
                    a: 1.5e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH2(S) + O <=> H + HCO",
                ArrheniusParams {
                    a: 1.5e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH3 + O <=> CH2O + H",
                ArrheniusParams {
                    a: 5.06e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH4 + O <=> CH3 + OH",
                ArrheniusParams {
                    a: 1.02e9,
                    n: 1.5,
                    ea: 8600.0,
                },
            ),
            (
                "HCO + O <=> CO + OH",
                ArrheniusParams {
                    a: 3.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "HCO + O <=> CO2 + H",
                ArrheniusParams {
                    a: 3.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH2O + O <=> HCO + OH",
                ArrheniusParams {
                    a: 3.9e13,
                    n: 0.0,
                    ea: 3540.0,
                },
            ),
            (
                "CH2OH + O <=> CH2O + OH",
                ArrheniusParams {
                    a: 1.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH3O + O <=> CH2O + OH",
                ArrheniusParams {
                    a: 1.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH3OH + O <=> CH3O + OH",
                ArrheniusParams {
                    a: 130000.0,
                    n: 2.5,
                    ea: 20920.0,
                },
            ),
            (
                "C2H + O <=> CH + CO",
                ArrheniusParams {
                    a: 5.0e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "C2H2 + O <=> H + HCCO",
                ArrheniusParams {
                    a: 1.35e7,
                    n: 2.0,
                    ea: 1900.0,
                },
            ),
            (
                "C2H2 + O <=> CH2 + CO",
                ArrheniusParams {
                    a: 6.94e6,
                    n: 2.0,
                    ea: 1900.0,
                },
            ),
            (
                "C2H3 + O <=> H + CH2CO",
                ArrheniusParams {
                    a: 4.8e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "C2H4 + O <=> CH3 + HCO",
                ArrheniusParams {
                    a: 1.25e7,
                    n: 1.83,
                    ea: 220.0,
                },
            ),
            (
                "C2H5 + O <=> CH2O + CH3",
                ArrheniusParams {
                    a: 2.24e13,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "C2H6 + O <=> C2H5 + OH",
                ArrheniusParams {
                    a: 8.98e7,
                    n: 1.92,
                    ea: 5690.0,
                },
            ),
            (
                "HCCO + O <=> H + 2 CO",
                ArrheniusParams {
                    a: 1.0e14,
                    n: 0.0,
                    ea: 0.0,
                },
            ),
            (
                "CH2CO + O <=> HCCO + OH",
                ArrheniusParams {
                    a: 1.0e13,
                    n: 0.0,
                    ea: 8000.0,
                },
            ),
            (
                "CH2CO + O <=> CH2 + CO2",
                ArrheniusParams {
                    a: 1.75e12,
                    n: 0.0,
                    ea: 1350.0,
                },
            ),
        ];

        for (equation, arrhenius) in sample_reactions {
            let reaction = Reaction {
                equation: equation.to_string(),
                coefficients: vec![1.0, 1.0, 1.0, 1.0],
                arrhenius_params: arrhenius,
            };
            self.reactions.insert(equation.to_string(), reaction);
            self.reaction_keys.push(equation.to_string());
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Выбор и добавление реакций и механизмов")
            .open(open)
            .default_size([1200.0, 800.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    // Left panel - Reaction list
                    ui.vertical(|ui| {
                        ui.set_width(400.0);
                        ui.heading("Список реакций");
                        // Search filter
                        ui.horizontal(|ui| {
                            ui.label("Поиск:");
                            ui.text_edit_singleline(&mut self.search_filter);
                        });
                        ui.separator();
                        egui::ScrollArea::vertical()
                            .max_height(500.0)
                            .show(ui, |ui| {
                                for reaction_key in &self.reaction_keys {
                                    if let Some(reaction) = self.reactions.get(reaction_key) {
                                        if self.search_filter.is_empty() ||
                                           reaction.equation.to_lowercase().contains(&self.search_filter.to_lowercase()) {
                                            let is_selected = self.selected_reaction_key.as_ref() == Some(reaction_key);
                                            if ui.selectable_label(is_selected, &reaction.equation).clicked() {
                                                self.selected_reaction_key = Some(reaction_key.clone());
                                            }
                                        }
                                    }
                                }
                            });
                        ui.separator();
                        // Dropdown for mechanism selection
                        egui::ComboBox::from_label("Механизм")
                            .selected_text("Cantera")
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut (), (), "Cantera");
                                ui.selectable_value(&mut (), (), "CHEMKIN");
                                ui.selectable_value(&mut (), (), "Custom");
                            });
                    });
                    ui.separator();
                    // Right panel - Reaction details and controls
                    ui.vertical(|ui| {
                        ui.heading("Данные выбранной реакции");
                        if let Some(reaction_key) = &self.selected_reaction_key {
                            if let Some(reaction) = self.reactions.get(reaction_key) {
                                // Display selected reaction details
                                ui.group(|ui| {
                                    ui.set_min_height(150.0);
                                    ui.label(format!("Реакция: {}", reaction.equation));
                                    ui.label("Параметры Аррениуса:");
                                    ui.label(format!("A: {:.2e}", reaction.arrhenius_params.a));
                                    ui.label(format!("n: {:.1}", reaction.arrhenius_params.n));
                                    ui.label(format!("Ea: {:.1}", reaction.arrhenius_params.ea));
                                    ui.separator();
                                    ui.label("Коэффициенты:");
                                    for (i, coeff) in reaction.coefficients.iter().enumerate() {
                                        ui.label(format!("k{}: {:.6}", i + 1, coeff));
                                    }
                                });
                            }
                        } else {
                            ui.group(|ui| {
                                ui.set_min_height(150.0);
                                ui.label("Выберите реакцию из списка");
                            });
                        }
                        ui.separator();
                        // Action buttons
                        ui.horizontal(|ui| {
                            if ui.button("Записать реакцию для расчета").clicked() {
                                println!("Saving reaction for calculation");
                            }
                            if ui.button("взять все реакции из механизма").clicked() {
                                println!("Taking all reactions from mechanism");
                            }
                        });
                        ui.horizontal(|ui| {
                            if ui.button("Искать").clicked() {
                                println!("Searching reactions");
                            }
                            if ui.button("построить подмеханизм").clicked() {
                                println!("Building sub-mechanism");
                            }
                        });
                        ui.horizontal(|ui| {
                            if ui.button("подмеханизм в расчет").clicked() {
                                println!("Adding sub-mechanism to calculation");
                            }
                        });
                        ui.separator();
                        // Input section
                        ui.heading("Построить подмеханизм для данных реагентов");                       
                        ui.horizontal(|ui| {
                            ui.label("Введите вещество для поиска:");
                            ui.text_edit_singleline(&mut self.mechanism_input);
                        });
                        // Radio buttons for reaction type
                        ui.horizontal(|ui| {
                            ui.radio_value(&mut self.reaction_type, ReactionType::Mechanism, "Механизм");
                            ui.radio_value(&mut self.reaction_type, ReactionType::Reaction, "Реакция");
                        });
                        ui.horizontal(|ui| {
                            ui.label("Введите новую реакцию или механизм:");
                            ui.text_edit_singleline(&mut self.reaction_input);
                            if ui.button("Записать").clicked() {
                                if !self.reaction_input.is_empty() {
                                    let new_reaction = Reaction {
                                        equation: self.reaction_input.clone(),
                                        coefficients: vec![1.0, 1.0, 1.0, 1.0],
                                        arrhenius_params: ArrheniusParams {
                                            a: 1.0e13,
                                            n: 0.0,
                                            ea: 0.0,
                                        },
                                    };
                                    self.reactions.insert(self.reaction_input.clone(), new_reaction);
                                    self.reaction_keys.push(self.reaction_input.clone());
                                    self.reaction_input.clear();
                                }
                            }
                        });
                        ui.separator();
                        // Bottom section
                        ui.label("Выберите из выпадающего списка в левой нижней углу механизм в который желаете");
                        ui.label("записать реакцию. Реакция записывается в формате 1:{eq: CH3+C2H4 <=> C2H3+CH4,");
                        ui.label("'type': 'elem', 'react': {CH3: 1, 'C2H4': 1}, 'Arrhenius': [227000.0, 2.0, 9200.0]}");
                        if ui.button("Записать из файла").clicked() {
                            println!("Importing from file");
                        }
                    });
                });
            });
    }
}
