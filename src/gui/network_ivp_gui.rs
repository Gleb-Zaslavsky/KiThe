use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::time::Duration;

use egui::{Color32, ComboBox, RichText};
use nalgebra::{DMatrix, DVector};
use strum::IntoEnumIterator;

use RustedSciThe::numerical::StiffAPI::{
    StiffCallbackPolicy, StiffFrontend, StiffMatrixLayout, StiffSolverKind,
};
use RustedSciThe::numerical::logging::SolverLogLevel;
use RustedSciThe::numerical::{CancellationToken, SolverTelemetryMode};
use RustedSciThe::symbolic::symbolic_engine::Expr;

use crate::Kinetics::reaction_network_ivp::{
    KineticOrder, NetworkIvpOptions, NetworkIvpSolution, NetworkProblem, RateCoefficient, RateLaw,
    ReactionSpec, SpeciesSpec, StoichTerm, TemperatureBinding, TemperatureProgram,
};
use crate::Kinetics::solid_state_kinetics_IVP::{KineticModel, KineticModelNames};
use crate::gui::gui_plot::PlotWindow;

#[derive(Clone, Debug)]
struct SpeciesDraft {
    id: String,
    initial: String,
}

#[derive(Clone, Debug)]
struct TermDraft {
    species: String,
    coefficient: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RateKind {
    MassAction,
    Empirical,
    Custom,
    SolidState,
}

#[derive(Clone, Debug)]
struct ReactionDraft {
    id: String,
    reactants: Vec<TermDraft>,
    products: Vec<TermDraft>,
    kind: RateKind,
    arrhenius: bool,
    prefactor: String,
    energy: String,
    orders: Vec<TermDraft>,
    custom: String,
    model: KineticModelNames,
    model_params: Vec<String>,
    capacity: String,
    initial_conversion: String,
}

impl ReactionDraft {
    fn new(id: String, reactant: String, product: String) -> Self {
        Self {
            id,
            reactants: vec![TermDraft {
                species: reactant,
                coefficient: "1".into(),
            }],
            products: vec![TermDraft {
                species: product,
                coefficient: "1".into(),
            }],
            kind: RateKind::MassAction,
            arrhenius: false,
            prefactor: "0.5".into(),
            energy: "5000".into(),
            orders: Vec::new(),
            custom: "0.5*c0".into(),
            model: KineticModelNames::F1,
            model_params: Vec::new(),
            capacity: "1".into(),
            initial_conversion: "0.1".into(),
        }
    }

    fn build(&self) -> Result<ReactionSpec, String> {
        let coefficient = || {
            Ok::<_, String>(if self.arrhenius {
                RateCoefficient::Arrhenius {
                    prefactor: number(&self.prefactor, "prefactor")?,
                    activation_energy: number(&self.energy, "activation energy")?,
                }
            } else {
                RateCoefficient::Constant(number(&self.prefactor, "rate coefficient")?)
            })
        };
        let rate_law = match self.kind {
            RateKind::MassAction => RateLaw::MassAction {
                coefficient: coefficient()?,
            },
            RateKind::Empirical => RateLaw::EmpiricalMassAction {
                coefficient: coefficient()?,
                orders: self
                    .orders
                    .iter()
                    .map(|term| {
                        Ok(KineticOrder::new(
                            term.species.trim(),
                            number(&term.coefficient, "kinetic order")?,
                        ))
                    })
                    .collect::<Result<Vec<_>, String>>()?,
            },
            RateKind::Custom => RateLaw::CustomExpr(
                Expr::try_parse_expression(self.custom.trim())
                    .map_err(|error| format!("reaction '{}': {error}", self.id))?,
            ),
            RateKind::SolidState => {
                let names = self.model.required_params();
                let params = names
                    .iter()
                    .enumerate()
                    .map(|(i, name)| {
                        number(self.model_params.get(i).map_or("", String::as_str), name)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let model = KineticModel::from_name_and_params(self.model, params)
                    .map_err(|error| format!("reaction '{}': {error}", self.id))?;
                RateLaw::SolidState {
                    coefficient: coefficient()?,
                    model,
                    capacity: number(&self.capacity, "capacity")?,
                    initial_conversion: number(&self.initial_conversion, "initial conversion")?,
                }
            }
        };
        Ok(ReactionSpec {
            id: self.id.trim().into(),
            reactants: build_terms(&self.reactants)?,
            products: build_terms(&self.products)?,
            rate_law,
        })
    }
}

fn number(value: &str, field: &str) -> Result<f64, String> {
    value
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("{field}: expected a number"))
}

fn build_terms(terms: &[TermDraft]) -> Result<Vec<StoichTerm>, String> {
    terms
        .iter()
        .map(|term| {
            Ok(StoichTerm::new(
                term.species.trim(),
                number(&term.coefficient, "stoichiometric coefficient")?,
            ))
        })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Page {
    Physics,
    Solver,
    Result,
    Help,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HelpLanguage {
    Russian,
    English,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TemperatureKind {
    Constant,
    Linear,
    Symbolic,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum LayoutKind {
    Dense,
    Sparse,
    BandedAuto,
    Banded,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum CallbackKind {
    Default,
    Sequential,
    Auto,
    Parallel,
}

struct Running {
    receiver: Receiver<Result<Vec<NetworkIvpSolution>, String>>,
    cancellation: CancellationToken,
}

pub struct NetworkIvpApp {
    page: Page,
    help_language: HelpLanguage,
    species: Vec<SpeciesDraft>,
    reactions: Vec<ReactionDraft>,
    temperature_kind: TemperatureKind,
    temperature: String,
    heating_rate: String,
    symbolic_temperature: String,
    start_time: String,
    end_time: String,
    family: StiffSolverKind,
    frontend: StiffFrontend,
    layout: LayoutKind,
    lower_band: String,
    upper_band: String,
    callback: CallbackKind,
    min_work: String,
    rtol: String,
    atol: String,
    max_step: String,
    telemetry: SolverTelemetryMode,
    logging: SolverLogLevel,
    running: Option<Running>,
    result: Option<NetworkIvpSolution>,
    plot: Option<PlotWindow>,
    message: Option<String>,
    series_mode: bool,
    series_values: String,
    series_id_prefix: String,
    completed_series: Option<Vec<NetworkIvpSolution>>,
    publish_requested: bool,
    last_run_description: String,
    existing_experiment_ids: std::collections::HashSet<String>,
}

impl Default for NetworkIvpApp {
    fn default() -> Self {
        Self {
            page: Page::Physics,
            help_language: HelpLanguage::Russian,
            species: vec![
                SpeciesDraft {
                    id: "A".into(),
                    initial: "1".into(),
                },
                SpeciesDraft {
                    id: "B".into(),
                    initial: "0".into(),
                },
            ],
            reactions: vec![ReactionDraft::new("r1".into(), "A".into(), "B".into())],
            temperature_kind: TemperatureKind::Constant,
            temperature: "700".into(),
            heating_rate: "1".into(),
            symbolic_temperature: "700+t".into(),
            start_time: "0".into(),
            end_time: "2".into(),
            family: StiffSolverKind::Bdf,
            frontend: StiffFrontend::ExprLegacy,
            layout: LayoutKind::Dense,
            lower_band: "1".into(),
            upper_band: "1".into(),
            callback: CallbackKind::Default,
            min_work: "64".into(),
            rtol: "1e-8".into(),
            atol: "1e-10".into(),
            max_step: "0.1".into(),
            telemetry: SolverTelemetryMode::Off,
            logging: SolverLogLevel::Off,
            running: None,
            result: None,
            plot: None,
            message: None,
            series_mode: false,
            series_values: "400, 500, 600".into(),
            series_id_prefix: "network".into(),
            completed_series: None,
            publish_requested: false,
            last_run_description: String::new(),
            existing_experiment_ids: std::collections::HashSet::new(),
        }
    }
}

impl NetworkIvpApp {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_for_experimental_kinetics() -> Self {
        Self {
            series_mode: true,
            ..Self::default()
        }
    }

    pub fn take_publish_request(&mut self) -> Option<(String, String, Vec<NetworkIvpSolution>)> {
        if !std::mem::take(&mut self.publish_requested) {
            return None;
        }
        self.completed_series.take().map(|solutions| {
            (
                self.series_id_prefix.trim().to_string(),
                self.last_run_description.clone(),
                solutions,
            )
        })
    }

    pub fn completed_run_count(&self) -> Option<usize> {
        self.completed_series.as_ref().map(Vec::len)
    }

    pub fn set_existing_experiment_ids(&mut self, ids: impl IntoIterator<Item = String>) {
        self.existing_experiment_ids = ids.into_iter().collect();
    }

    fn build_problem(&self) -> Result<NetworkProblem, String> {
        let species = self
            .species
            .iter()
            .map(|entry| {
                Ok(SpeciesSpec {
                    id: entry.id.trim().into(),
                    initial_concentration: number(&entry.initial, "initial concentration")?,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let reactions = self
            .reactions
            .iter()
            .map(ReactionDraft::build)
            .collect::<Result<Vec<_>, _>>()?;
        let temperature = match self.temperature_kind {
            TemperatureKind::Constant => {
                TemperatureProgram::Constant(number(&self.temperature, "temperature")?)
            }
            TemperatureKind::Linear => TemperatureProgram::Linear {
                initial_kelvin: number(&self.temperature, "temperature at t=0")?,
                kelvin_per_second: number(&self.heating_rate, "heating rate")?,
            },
            TemperatureKind::Symbolic => TemperatureProgram::Symbolic(
                Expr::try_parse_expression(self.symbolic_temperature.trim())
                    .map_err(|error| format!("temperature expression: {error}"))?,
            ),
        };
        NetworkProblem::new(
            species,
            reactions,
            temperature,
            number(&self.start_time, "start time")?,
            number(&self.end_time, "end time")?,
        )
        .map_err(|error| error.to_string())
    }

    fn build_options(&self) -> Result<NetworkIvpOptions, String> {
        let matrix_layout = match self.layout {
            LayoutKind::Dense => StiffMatrixLayout::Dense,
            LayoutKind::Sparse => StiffMatrixLayout::Sparse,
            LayoutKind::BandedAuto => StiffMatrixLayout::BandedAuto,
            LayoutKind::Banded => StiffMatrixLayout::Banded {
                lower: self
                    .lower_band
                    .trim()
                    .parse()
                    .map_err(|_| "lower bandwidth: expected integer")?,
                upper: self
                    .upper_band
                    .trim()
                    .parse()
                    .map_err(|_| "upper bandwidth: expected integer")?,
            },
        };
        let callback_policy = match if self.family == StiffSolverKind::BackwardEuler {
            CallbackKind::Default
        } else {
            self.callback
        } {
            CallbackKind::Default => None,
            CallbackKind::Sequential => Some(StiffCallbackPolicy::Sequential),
            CallbackKind::Auto | CallbackKind::Parallel => {
                let min_work = self
                    .min_work
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| "callback threshold: expected integer")?;
                Some(if self.callback == CallbackKind::Auto {
                    StiffCallbackPolicy::Auto { min_work }
                } else {
                    StiffCallbackPolicy::Parallel { min_work }
                })
            }
        };
        Ok(NetworkIvpOptions {
            family: self.family,
            frontend: self.frontend,
            matrix_layout,
            callback_policy,
            rtol: if self.family == StiffSolverKind::BackwardEuler {
                1e-8
            } else {
                number(&self.rtol, "relative tolerance")?
            },
            atol: if self.family == StiffSolverKind::BackwardEuler {
                1e-10
            } else {
                number(&self.atol, "absolute tolerance")?
            },
            max_step: number(&self.max_step, "maximum step")?,
            telemetry: self.telemetry,
            logging: self.logging,
        })
    }

    fn start(&mut self) {
        let request = self
            .build_problem()
            .and_then(|problem| self.build_options().map(|options| (problem, options)));
        let (problem, options) = match request {
            Ok(request) => request,
            Err(error) => {
                self.message = Some(error);
                return;
            }
        };
        self.last_run_description = format!(
            "species={:?}; reactions={:?}; temperature_program={:?}; interval=[{}, {}] s",
            self.species, self.reactions, self.temperature_kind, self.start_time, self.end_time
        );
        if self.series_mode {
            let bindings = match self.parse_series_bindings() {
                Ok(bindings) => bindings,
                Err(error) => {
                    self.message = Some(error);
                    return;
                }
            };
            let initial = problem.initial_concentrations();
            let (sender, receiver) = mpsc::channel();
            std::thread::spawn(move || {
                let result = (|| {
                    let mut prepared =
                        crate::Kinetics::reaction_network_ivp::PreparedNetworkIvp::new(
                            &problem, &options,
                        )
                        .map_err(|error| error.to_string())?;
                    let mut solutions = Vec::with_capacity(bindings.len());
                    for binding in bindings {
                        let solution = prepared
                            .solve_with_temperature_binding(&initial, binding)
                            .map(Clone::clone)
                            .map_err(|error| error.to_string())?;
                        solutions.push(solution);
                    }
                    Ok(solutions)
                })();
                let _ = sender.send(result);
            });
            self.completed_series = None;
            self.running = Some(Running {
                receiver,
                cancellation: CancellationToken::new(),
            });
            self.message = Some(format!(
                "Running {} conditions",
                self.parse_series_bindings().map_or(0, |v| v.len())
            ));
            return;
        }
        let (sender, receiver) = mpsc::channel();
        let cancellation = CancellationToken::new();
        let worker_token = cancellation.clone();
        std::thread::spawn(move || {
            let result = (|| {
                let mut prepared = crate::Kinetics::reaction_network_ivp::PreparedNetworkIvp::new(
                    &problem, &options,
                )
                .map_err(|error| error.to_string())?;
                let initial = problem.initial_concentrations();
                let solution = prepared
                    .solve_with_cancellation_token(&initial, worker_token)
                    .map(Clone::clone)
                    .map_err(|error| error.to_string())?;
                Ok(vec![solution])
            })();
            let _ = sender.send(result);
        });
        self.running = Some(Running {
            receiver,
            cancellation,
        });
        self.message = Some("Running".into());
    }

    fn poll(&mut self, ctx: &egui::Context) {
        let Some(running) = &self.running else {
            return;
        };
        match running.receiver.try_recv() {
            Ok(Ok(solutions)) => {
                self.message = Some("Completed".into());
                self.result = solutions.last().cloned();
                self.completed_series = self.series_mode.then_some(solutions);
                self.running = None;
                self.page = Page::Result;
            }
            Ok(Err(error)) => {
                self.message = Some(error);
                self.running = None;
            }
            Err(TryRecvError::Disconnected) => {
                self.message = Some("Solver worker stopped unexpectedly".into());
                self.running = None;
            }
            Err(TryRecvError::Empty) => ctx.request_repaint_after(Duration::from_millis(80)),
        }
    }

    fn parse_series_bindings(&self) -> Result<Vec<TemperatureBinding>, String> {
        let values = parse_finite_csv(&self.series_values)?;
        if values.is_empty() {
            return Err("Enter at least one temperature or heating rate".into());
        }
        if values.len() > 100 {
            return Err("A batch is limited to 100 conditions".into());
        }
        let mut unique = std::collections::HashSet::with_capacity(values.len());
        if values.iter().any(|value| !unique.insert(value.to_bits())) {
            return Err("Series conditions must not contain duplicates".into());
        }
        let bindings: Vec<TemperatureBinding> = match self.temperature_kind {
            TemperatureKind::Constant => {
                if values.iter().any(|value| *value <= 0.0) {
                    return Err("Temperatures must be greater than 0 K".into());
                }
                values
                    .into_iter()
                    .map(TemperatureBinding::Constant)
                    .collect()
            }
            TemperatureKind::Linear => {
                let initial_kelvin = number(&self.temperature, "temperature at t=0")?;
                if initial_kelvin <= 0.0 || values.iter().any(|rate| *rate <= 0.0) {
                    return Err("Initial temperature and heating rates must be positive".into());
                }
                values
                    .into_iter()
                    .map(|rate| TemperatureBinding::Linear {
                        initial_kelvin,
                        kelvin_per_second: rate / 60.0,
                    })
                    .collect()
            }
            TemperatureKind::Symbolic => {
                return Err(
                    "Condition series supports constant temperature or linear heating only".into(),
                );
            }
        };
        let mut ids = std::collections::HashSet::with_capacity(bindings.len());
        for binding in &bindings {
            let id = network_condition_experiment_id(&self.series_id_prefix, *binding);
            if self.existing_experiment_ids.contains(&id) || !ids.insert(id.clone()) {
                return Err(format!(
                    "Experiment ID '{id}' already exists or repeats in this batch"
                ));
            }
        }
        Ok(bindings)
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        self.poll(ctx);
        egui::Window::new("Reaction network IVP")
            .open(open)
            .default_size([950.0, 680.0])
            .min_size([620.0, 420.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut self.page, Page::Physics, "Task");
                    ui.selectable_value(&mut self.page, Page::Solver, "Solver");
                    ui.selectable_value(&mut self.page, Page::Result, "Result");
                    ui.selectable_value(&mut self.page, Page::Help, "Help");
                });
                ui.separator();
                egui::ScrollArea::vertical()
                    .max_height(ui.available_height() - 55.0)
                    .show(ui, |ui| match self.page {
                        Page::Physics => self.physics_ui(ui),
                        Page::Solver => self.solver_ui(ui),
                        Page::Result => self.result_ui(ui),
                        Page::Help => self.help_ui(ui),
                    });
                ui.separator();
                ui.horizontal(|ui| {
                    if self.running.is_none() {
                        if ui.button("Run").clicked() {
                            self.start();
                        }
                    } else if self.series_mode {
                        ui.add_enabled(false, egui::Button::new("Running condition series..."));
                    } else if ui.button("Cancel").clicked() {
                        if let Some(running) = &self.running {
                            running.cancellation.cancel();
                        }
                    }
                    if let Some(message) = &self.message {
                        let color = if message == "Completed" || message == "Running" {
                            ui.visuals().text_color()
                        } else {
                            Color32::LIGHT_RED
                        };
                        ui.colored_label(color, message);
                    }
                });
            });
        if !*open {
            if let Some(running) = &self.running {
                running.cancellation.cancel();
            }
        }
        if let Some(plot) = &mut self.plot {
            plot.show(ctx);
        }
    }

    fn physics_ui(&mut self, ui: &mut egui::Ui) {
        ui.heading("Conditions");
        ui.add_enabled(
            self.running.is_none(),
            egui::Checkbox::new(&mut self.series_mode, "Calculate a condition series"),
        );
        ui.horizontal(|ui| {
            labeled_edit(ui, "From, s", &mut self.start_time, 90.0);
            labeled_edit(ui, "To, s", &mut self.end_time, 90.0);
        });
        if self.series_mode {
            ui.horizontal(|ui| {
                let label = ui.label(match self.temperature_kind {
                    TemperatureKind::Constant => "Temperatures, K:",
                    TemperatureKind::Linear => "Heating rates, K/min:",
                    TemperatureKind::Symbolic => "Conditions:",
                });
                ui.text_edit_singleline(&mut self.series_values)
                    .labelled_by(label.id);
                let label = ui.label("Experiment ID prefix:");
                ui.text_edit_singleline(&mut self.series_id_prefix)
                    .labelled_by(label.id);
            });
            ui.label("Each condition starts from the same initial concentrations and uses the shared time interval.");
            match self.parse_series_bindings() {
                Ok(bindings) => {
                    ui.label(format!("{} experiments will be created:", bindings.len()));
                    for binding in bindings {
                        ui.label(network_condition_experiment_id(
                            &self.series_id_prefix,
                            binding,
                        ));
                    }
                }
                Err(error) => {
                    ui.colored_label(Color32::LIGHT_RED, error);
                }
            }
        }
        ui.horizontal(|ui| {
            let label = ui.label("Temperature");
            ComboBox::from_id_salt("network_temperature_kind")
                .selected_text(match self.temperature_kind {
                    TemperatureKind::Constant => "Constant",
                    TemperatureKind::Linear => "Linear",
                    TemperatureKind::Symbolic => "T(t)",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.temperature_kind,
                        TemperatureKind::Constant,
                        "Constant",
                    );
                    ui.selectable_value(
                        &mut self.temperature_kind,
                        TemperatureKind::Linear,
                        "Linear",
                    );
                    ui.selectable_value(
                        &mut self.temperature_kind,
                        TemperatureKind::Symbolic,
                        "T(t)",
                    );
                })
                .response
                .labelled_by(label.id);
            match self.temperature_kind {
                TemperatureKind::Constant => {
                    if !self.series_mode {
                        labeled_edit(ui, "K", &mut self.temperature, 100.0);
                    }
                }
                TemperatureKind::Linear => {
                    labeled_edit(ui, "T(0), K", &mut self.temperature, 90.0);
                    if !self.series_mode {
                        labeled_edit(ui, "K/s", &mut self.heating_rate, 90.0);
                    }
                }
                TemperatureKind::Symbolic => {
                    let label = ui.label("T(t), K");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.symbolic_temperature)
                            .desired_width(250.0),
                    )
                    .labelled_by(label.id);
                }
            }
        });
        ui.add_space(10.0);
        ui.horizontal(|ui| {
            ui.heading("Species");
            if ui.button("Add species").clicked() {
                self.species.push(SpeciesDraft {
                    id: format!("S{}", self.species.len() + 1),
                    initial: "0".into(),
                });
            }
        });
        let mut remove_species = None;
        egui::Grid::new("network_species")
            .striped(true)
            .show(ui, |ui| {
                ui.strong("ID");
                ui.strong("Initial, mol/m^3");
                ui.end_row();
                for (index, entry) in self.species.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        labeled_edit(
                            ui,
                            &format!("Species {} ID", index + 1),
                            &mut entry.id,
                            150.0,
                        );
                    });
                    ui.horizontal(|ui| {
                        labeled_edit(ui, "Initial", &mut entry.initial, 130.0);
                    });
                    if ui.button("-").on_hover_text("Remove species").clicked() {
                        remove_species = Some(index);
                    }
                    ui.end_row();
                }
            });
        if let Some(index) = remove_species {
            self.species.remove(index);
        }
        ui.add_space(12.0);
        ui.horizontal(|ui| {
            ui.heading("Reactions");
            if ui.button("Add reaction").clicked() {
                let from = self.species.first().map_or("A", |s| &s.id).to_owned();
                let to = self.species.get(1).map_or("B", |s| &s.id).to_owned();
                self.reactions.push(ReactionDraft::new(
                    format!("r{}", self.reactions.len() + 1),
                    from,
                    to,
                ));
            }
        });
        let mut remove_reaction = None;
        for (index, reaction) in self.reactions.iter_mut().enumerate() {
            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.horizontal(|ui| {
                    labeled_edit(ui, "ID", &mut reaction.id, 130.0);
                    if ui.button("-").on_hover_text("Remove reaction").clicked() {
                        remove_reaction = Some(index);
                    }
                });
                terms_ui(
                    ui,
                    &mut reaction.reactants,
                    "Reactants",
                    &self.species,
                    index,
                );
                terms_ui(ui, &mut reaction.products, "Products", &self.species, index);
                ui.horizontal(|ui| {
                    let label = ui.label("Rate law");
                    ComboBox::from_id_salt(("network_law", index))
                        .selected_text(match reaction.kind {
                            RateKind::MassAction => "Mass action",
                            RateKind::Empirical => "Empirical orders",
                            RateKind::Custom => "Expression",
                            RateKind::SolidState => "Solid state",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut reaction.kind,
                                RateKind::MassAction,
                                "Mass action",
                            );
                            ui.selectable_value(
                                &mut reaction.kind,
                                RateKind::Empirical,
                                "Empirical orders",
                            );
                            ui.selectable_value(&mut reaction.kind, RateKind::Custom, "Expression");
                            ui.selectable_value(
                                &mut reaction.kind,
                                RateKind::SolidState,
                                "Solid state",
                            );
                        })
                        .response
                        .labelled_by(label.id);
                });
                match reaction.kind {
                    RateKind::Custom => {
                        ui.horizontal(|ui| {
                            let label = ui.label("r(c,T,t)");
                            ui.add(
                                egui::TextEdit::singleline(&mut reaction.custom)
                                    .desired_width(300.0),
                            )
                            .labelled_by(label.id)
                            .on_hover_text("Use c0, c1, ... in species order, T and t");
                        });
                    }
                    _ => {
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut reaction.arrhenius, "Arrhenius");
                            labeled_edit(
                                ui,
                                if reaction.arrhenius { "A" } else { "k" },
                                &mut reaction.prefactor,
                                100.0,
                            );
                            if reaction.arrhenius {
                                labeled_edit(ui, "Ea, J/mol", &mut reaction.energy, 100.0);
                            }
                        });
                        if reaction.kind == RateKind::Empirical {
                            terms_ui(
                                ui,
                                &mut reaction.orders,
                                "Kinetic orders",
                                &self.species,
                                index,
                            );
                        }
                        if reaction.kind == RateKind::SolidState {
                            ui.horizontal(|ui| {
                                let label = ui.label("Model");
                                let before = reaction.model;
                                ComboBox::from_id_salt(("network_model", index))
                                    .selected_text(format!("{:?}", reaction.model))
                                    .show_ui(ui, |ui| {
                                        for model in KineticModelNames::iter() {
                                            ui.selectable_value(
                                                &mut reaction.model,
                                                model,
                                                format!("{model:?}"),
                                            );
                                        }
                                    })
                                    .response
                                    .labelled_by(label.id);
                                if reaction.model != before {
                                    reaction.model_params.clear();
                                }
                                labeled_edit(ui, "q, mol/m^3", &mut reaction.capacity, 85.0);
                                labeled_edit(
                                    ui,
                                    "alpha(0)",
                                    &mut reaction.initial_conversion,
                                    70.0,
                                );
                            });
                            let names = reaction.model.required_params();
                            reaction.model_params.resize(names.len(), "1".into());
                            if !names.is_empty() {
                                ui.horizontal(|ui| {
                                    for (name, value) in
                                        names.iter().zip(reaction.model_params.iter_mut())
                                    {
                                        labeled_edit(ui, name, value, 65.0);
                                    }
                                });
                            }
                        }
                    }
                }
            });
            ui.add_space(6.0);
        }
        if let Some(index) = remove_reaction {
            self.reactions.remove(index);
        }
    }

    fn solver_ui(&mut self, ui: &mut egui::Ui) {
        ui.heading("Numerical solver");
        egui::Grid::new("network_solver_settings")
            .num_columns(2)
            .spacing([16.0, 10.0])
            .show(ui, |ui| {
                let label = ui.label("Family");
                ComboBox::from_id_salt("network_family")
                    .selected_text(format!("{:?}", self.family))
                    .show_ui(ui, |ui| {
                        for (family, label) in [
                            (StiffSolverKind::Bdf, "BDF"),
                            (StiffSolverKind::Radau, "Radau"),
                            (StiffSolverKind::Lsode2, "LSODE2"),
                            (StiffSolverKind::BackwardEuler, "Backward Euler"),
                        ] {
                            ui.selectable_value(&mut self.family, family, label);
                        }
                    })
                    .response
                    .labelled_by(label.id);
                ui.end_row();
                let label = ui.label("Symbolic frontend");
                ComboBox::from_id_salt("network_frontend")
                    .selected_text(format!("{:?}", self.frontend))
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.frontend,
                            StiffFrontend::ExprLegacy,
                            "ExprLegacy",
                        );
                        ui.selectable_value(
                            &mut self.frontend,
                            StiffFrontend::AtomViewNative,
                            "AtomView",
                        );
                    })
                    .response
                    .labelled_by(label.id);
                ui.end_row();
                let label = ui.label("Matrix layout");
                ui.horizontal(|ui| {
                    ComboBox::from_id_salt("network_layout")
                        .selected_text(match self.layout {
                            LayoutKind::Dense => "Dense",
                            LayoutKind::Sparse => "Sparse",
                            LayoutKind::BandedAuto => "Banded auto",
                            LayoutKind::Banded => "Banded",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.layout, LayoutKind::Dense, "Dense");
                            ui.selectable_value(&mut self.layout, LayoutKind::Sparse, "Sparse");
                            ui.selectable_value(
                                &mut self.layout,
                                LayoutKind::BandedAuto,
                                "Banded auto",
                            );
                            ui.selectable_value(&mut self.layout, LayoutKind::Banded, "Banded");
                        })
                        .response
                        .labelled_by(label.id);
                    if self.layout == LayoutKind::Banded {
                        labeled_edit(ui, "lower", &mut self.lower_band, 45.0);
                        labeled_edit(ui, "upper", &mut self.upper_band, 45.0);
                    }
                });
                ui.end_row();
                let label = ui.label("Callbacks");
                ui.horizontal(|ui| {
                    if self.family == StiffSolverKind::BackwardEuler {
                        ui.label("Unavailable");
                    }
                    ui.add_enabled_ui(self.family != StiffSolverKind::BackwardEuler, |ui| {
                        ComboBox::from_id_salt("network_callback")
                            .selected_text(match self.callback {
                                CallbackKind::Default => "Default",
                                CallbackKind::Sequential => "Sequential",
                                CallbackKind::Auto => "Auto",
                                CallbackKind::Parallel => "Parallel",
                            })
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.callback,
                                    CallbackKind::Default,
                                    "Default",
                                );
                                ui.selectable_value(
                                    &mut self.callback,
                                    CallbackKind::Sequential,
                                    "Sequential",
                                );
                                ui.selectable_value(&mut self.callback, CallbackKind::Auto, "Auto");
                                ui.selectable_value(
                                    &mut self.callback,
                                    CallbackKind::Parallel,
                                    "Parallel",
                                );
                            })
                            .response
                            .labelled_by(label.id);
                        if matches!(self.callback, CallbackKind::Auto | CallbackKind::Parallel) {
                            labeled_edit(ui, "min work", &mut self.min_work, 65.0);
                        }
                    });
                });
                ui.end_row();
                for (label, value) in [
                    ("Relative tolerance", &mut self.rtol),
                    ("Absolute tolerance", &mut self.atol),
                ] {
                    let name = ui.label(label);
                    ui.add_enabled(
                        self.family != StiffSolverKind::BackwardEuler,
                        egui::TextEdit::singleline(value).desired_width(110.0),
                    )
                    .labelled_by(name.id);
                    ui.end_row();
                }
                let label = ui.label("Maximum step, s");
                ui.add(egui::TextEdit::singleline(&mut self.max_step).desired_width(110.0))
                    .labelled_by(label.id);
                ui.end_row();
                let label = ui.label("Telemetry");
                ComboBox::from_id_salt("network_telemetry")
                    .selected_text(format!("{:?}", self.telemetry))
                    .show_ui(ui, |ui| {
                        for mode in [
                            SolverTelemetryMode::Off,
                            SolverTelemetryMode::Counters,
                            SolverTelemetryMode::Timings,
                            SolverTelemetryMode::Detailed,
                        ] {
                            ui.selectable_value(&mut self.telemetry, mode, format!("{mode:?}"));
                        }
                    })
                    .response
                    .labelled_by(label.id);
                ui.end_row();
                let label = ui.label("Logging");
                ComboBox::from_id_salt("network_logging")
                    .selected_text(format!("{:?}", self.logging))
                    .show_ui(ui, |ui| {
                        for level in [
                            SolverLogLevel::Off,
                            SolverLogLevel::Error,
                            SolverLogLevel::Warn,
                            SolverLogLevel::Info,
                            SolverLogLevel::Debug,
                        ] {
                            ui.selectable_value(&mut self.logging, level, format!("{level:?}"));
                        }
                    })
                    .response
                    .labelled_by(label.id);
                ui.end_row();
            });
    }

    fn result_ui(&mut self, ui: &mut egui::Ui) {
        if let Some(solutions) = &self.completed_series {
            ui.heading(format!("Completed condition series ({})", solutions.len()));
            egui::Grid::new("network_completed_conditions")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Run");
                    ui.strong("T(0), K");
                    ui.strong("T(end), K");
                    ui.strong("Samples");
                    ui.end_row();
                    for (index, run) in solutions.iter().enumerate() {
                        ui.label(format!("{}", index + 1));
                        ui.label(format!(
                            "{:.4}",
                            run.temperature_kelvin.first().unwrap_or(&f64::NAN)
                        ));
                        ui.label(format!(
                            "{:.4}",
                            run.temperature_kelvin.last().unwrap_or(&f64::NAN)
                        ));
                        ui.label(run.time_seconds.len().to_string());
                        ui.end_row();
                    }
                });
            if ui.button("Add all runs to Experimental Kinetics").clicked() {
                self.publish_requested = true;
            }
            ui.separator();
        }
        let Some(solution) = &self.result else {
            ui.label("No result");
            return;
        };
        let Some(final_state) = solution.concentrations.last() else {
            return;
        };
        ui.heading("Final state");
        ui.label(format!(
            "t = {:.6e} s    T = {:.6e} K    samples = {}",
            solution.time_seconds.last().unwrap_or(&0.0),
            solution.temperature_kelvin.last().unwrap_or(&0.0),
            solution.time_seconds.len()
        ));
        egui::Grid::new("network_final_state")
            .striped(true)
            .show(ui, |ui| {
                ui.strong("Species");
                ui.strong("mol/m^3");
                ui.end_row();
                for (name, value) in solution.species_ids.iter().zip(final_state) {
                    ui.label(name);
                    ui.label(format!("{value:.8e}"));
                    ui.end_row();
                }
            });
        if let Some(conversions) = solution.conversions.last() {
            for (name, value) in solution.solid_state_reaction_ids.iter().zip(conversions) {
                ui.label(format!("{name}: alpha = {value:.8e}"));
            }
        }
        if ui.button("Plot concentrations").clicked() {
            let flat: Vec<f64> = solution
                .concentrations
                .iter()
                .flat_map(|row| row.iter().copied())
                .collect();
            self.plot = Some(PlotWindow::new(
                "Time, s".into(),
                solution.species_ids.clone(),
                DVector::from_vec(solution.time_seconds.clone()),
                DMatrix::from_row_slice(
                    solution.time_seconds.len(),
                    solution.species_ids.len(),
                    &flat,
                ),
            ));
        }
        if let Some(telemetry) = &solution.telemetry {
            ui.collapsing("Solver diagnostics", |ui| {
                ui.monospace(format!("{telemetry:#?}"));
            });
        }
        if let Some(diagnostics) = &solution.stage_diagnostics {
            ui.collapsing("KiThe stages", |ui| {
                ui.monospace(format!("{diagnostics:#?}"));
            });
        }
    }

    fn help_ui(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.help_language, HelpLanguage::Russian, "Русский");
            ui.selectable_value(&mut self.help_language, HelpLanguage::English, "English");
        });
        ui.separator();
        let sections: [(&str, &str); 4] = match self.help_language {
            HelpLanguage::Russian => [
                (
                    "Физическая модель",
                    "Решается система dC/dt = F(C, T(t)) для концентраций веществ. Температурная программа задаётся как постоянная, линейная или выражением T(t). Термодинамическое равновесие не рассчитывается.",
                ),
                (
                    "Доступные солверы",
                    "Доступны BDF, Radau, LSODE2 и Backward Euler из StiffAPI. Это семейства интеграторов для начальной задачи; доступность отдельных сочетаний фронтенда и линейной алгебры проверяется при подготовке.",
                ),
                (
                    "Основные настройки",
                    "На вкладке Solver выбираются семейство, символьный frontend, layout матрицы, callback policy, допуски, максимальный шаг, telemetry и logging. Backward Euler не использует общие rtol/atol и callback policy.",
                ),
                (
                    "Типичный workflow",
                    "На Task задайте вещества и начальные концентрации, реакции и закон скорости, температурную программу и интервал времени. Затем настройте Solver и нажмите Run. Результаты и графики концентраций доступны на Result; расчёт можно отменить.",
                ),
            ],
            HelpLanguage::English => [
                (
                    "Physical model",
                    "The solver integrates dC/dt = F(C, T(t)) for species concentrations. Temperature can be constant, linear, or defined by an expression T(t). Thermodynamic equilibrium is not calculated.",
                ),
                (
                    "Available solvers",
                    "BDF, Radau, LSODE2, and Backward Euler are available through StiffAPI. These are initial-value integrators; supported frontend and linear-algebra combinations are validated during preparation.",
                ),
                (
                    "Main settings",
                    "The Solver tab selects the family, symbolic frontend, matrix layout, callback policy, tolerances, maximum step, telemetry, and logging. Backward Euler does not use the common rtol/atol or callback policy.",
                ),
                (
                    "Typical workflow",
                    "On Task, define species and initial concentrations, reactions and rate laws, temperature program, and time interval. Configure Solver and press Run. Results and concentration plots are available on Result; a running solve can be cancelled.",
                ),
            ],
        };
        for (heading, body) in sections {
            ui.heading(heading);
            ui.label(body);
            ui.add_space(8.0);
        }
    }
}

fn parse_finite_csv(input: &str) -> Result<Vec<f64>, String> {
    let values = input
        .split([',', ';', ' '])
        .filter(|part| !part.trim().is_empty())
        .map(|part| number(part, "series condition"))
        .collect::<Result<Vec<_>, _>>()?;
    if values.iter().any(|value| !value.is_finite()) {
        return Err("Series conditions must be finite numbers".into());
    }
    Ok(values)
}

pub fn network_condition_experiment_id(prefix: &str, binding: TemperatureBinding) -> String {
    let label = match binding {
        TemperatureBinding::Constant(kelvin) => format!("{kelvin}K"),
        TemperatureBinding::Linear {
            kelvin_per_second, ..
        } => {
            format!("{}Kpm", kelvin_per_second * 60.0)
        }
    };
    let prefix = if prefix.trim().is_empty() {
        "network"
    } else {
        prefix.trim()
    };
    format!("{prefix}_{}", label.replace('.', "_"))
}

fn terms_ui(
    ui: &mut egui::Ui,
    terms: &mut Vec<TermDraft>,
    title: &str,
    species: &[SpeciesDraft],
    reaction_index: usize,
) {
    ui.horizontal(|ui| {
        ui.label(RichText::new(title).strong());
        if ui
            .button("+")
            .on_hover_text(format!("Add {title}"))
            .clicked()
        {
            terms.push(TermDraft {
                species: species.first().map_or("", |entry| entry.id.as_str()).into(),
                coefficient: "1".into(),
            });
        }
    });
    let mut remove = None;
    for (index, term) in terms.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            let name = ui.label(format!("{title} {} species", index + 1));
            ui.add(egui::TextEdit::singleline(&mut term.species).desired_width(130.0))
                .labelled_by(name.id);
            ComboBox::from_id_salt(("network_species_choice", reaction_index, title, index))
                .selected_text("v")
                .show_ui(ui, |ui| {
                    for entry in species {
                        if ui
                            .selectable_label(term.species == entry.id, &entry.id)
                            .clicked()
                        {
                            term.species = entry.id.clone();
                        }
                    }
                })
                .response
                .labelled_by(name.id);
            labeled_edit(ui, "Coefficient", &mut term.coefficient, 80.0);
            if ui
                .button("-")
                .on_hover_text(format!("Remove {title}"))
                .clicked()
            {
                remove = Some(index);
            }
        });
    }
    if let Some(index) = remove {
        terms.remove(index);
    }
}

fn labeled_edit(ui: &mut egui::Ui, label: &str, value: &mut String, width: f32) {
    let name = ui.label(label);
    ui.add(egui::TextEdit::singleline(value).desired_width(width))
        .labelled_by(name.id);
}

#[cfg(test)]
#[path = "network_ivp_gui_tests.rs"]
mod tests;
