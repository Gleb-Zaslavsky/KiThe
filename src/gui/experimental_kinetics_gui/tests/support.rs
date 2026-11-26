use crate::Kinetics::experimental_kinetics::experiment_series_main::TGASeries;
use crate::Kinetics::experimental_kinetics::math::ndarray_statistics::{
    FilterQualityReport, filter_quality_from_arrays,
};
use crate::Kinetics::experimental_kinetics::testing_mod::{
    AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind, SimulatedTGA,
    VirtualTGA,
};
use crate::gui::experimental_kinetics_gui::controller_buttons_and_panels::{
    NewExperimentDialogState, QuickActionPanelState,
};
use crate::gui::experimental_kinetics_gui::controller_filters::Mathematics;
use crate::gui::experimental_kinetics_gui::controller_kinetics::DirectProblemDialogState;
use crate::gui::experimental_kinetics_gui::controller_methods::KineticMethodsWindowState;
use crate::gui::experimental_kinetics_gui::controller_table::ColumnManagerState;
use crate::gui::experimental_kinetics_gui::experimental_kinetics_gui_main::PlotApp;
use crate::gui::experimental_kinetics_gui::kitheplot_wrapper::KiThePlotWindowState;
use crate::gui::experimental_kinetics_gui::model::PlotModel;
use crate::gui::experimental_kinetics_gui::test_options::TestOptions;
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use ndarray::Array1;
use polars::prelude::col;
use std::path::Path;

pub struct StoryState {
    pub app: PlotApp,
    open: bool,
}

pub fn empty_app() -> PlotApp {
    PlotApp {
        model: PlotModel::new(),
        quick_actions_state: QuickActionPanelState::default(),
        new_experiment_dialog: NewExperimentDialogState::new(),
        mathematics: Mathematics::new(),
        test_options: TestOptions::new(),
        column_manager_state: ColumnManagerState::new(),
        kithe_plot_window: KiThePlotWindowState::default(),
        kinetic_methods_state: KineticMethodsWindowState::default(),
        direct_problem_state: DirectProblemDialogState::new(),
        keyboard_input_focused_last_frame: false,
    }
}

pub fn with_app(app: PlotApp, run: impl FnOnce(&mut Harness<StoryState>)) {
    let state = StoryState { app, open: true };
    let mut harness = Harness::new_ui_state(
        |ui, state: &mut StoryState| state.app.show(ui.ctx(), &mut state.open),
        state,
    );
    run(&mut harness);
}

pub fn click_menu_entry<State>(harness: &mut Harness<State>, menu: &str, entry: &str) {
    harness.get_by_label(menu).click();
    harness.run();
    harness.get_by_label(entry).click();
    harness.run();
}

pub fn deterministic_config(isothermal: bool, seed: u64) -> AdvancedTGAConfig {
    AdvancedTGAConfig {
        n_points: 600,
        dt: 2.0,
        experiment_mode: if isothermal {
            ExperimentMode::Isothermal {
                temperatures: vec![500.0, 530.0, 560.0],
            }
        } else {
            ExperimentMode::NonIsothermal {
                t0: 300.0,
                heating_rates: vec![2.0, 5.0, 10.0],
            }
        },
        kinetic_model: KineticModel::ArrheniusSingle {
            m0: 100.0,
            k0: 2.0e5,
            e: 85_000.0,
            r: 8.314,
        },
        mass_noise: Some(NoiseConfig {
            kind: NoiseKind::Gaussian { sigma: 0.15 },
        }),
        temp_noise: Some(NoiseConfig {
            kind: NoiseKind::Gaussian { sigma: 0.08 },
        }),
        spikes: None,
        seed,
    }
}

pub fn write_simulated_csv(path: &Path, simulated: &SimulatedTGA, voltage: bool) {
    let mut csv = String::from("t,T,m\n");
    for ((time, temperature), mass) in simulated
        .virtual_tga
        .time
        .iter()
        .zip(&simulated.virtual_tga.temperature)
        .zip(&simulated.virtual_tga.mass)
    {
        // The default calibration in test settings is mass = -voltage + 40.
        let mass_or_voltage = if voltage { 40.0 - mass } else { *mass };
        csv.push_str(&format!(
            "{time:.8},{temperature:.8},{mass_or_voltage:.8}\n"
        ));
    }
    std::fs::write(path, csv).expect("write synthetic CSV fixture");
}

pub fn generated(config: &AdvancedTGAConfig) -> Vec<SimulatedTGA> {
    VirtualTGA::generate_series(config).expect("valid synthetic story-test configuration")
}

pub fn import_csv(
    harness: &mut Harness<StoryState>,
    path: &Path,
    experiment_id: &str,
    isothermal_temperature: Option<f64>,
    heating_rate: Option<f64>,
) {
    harness.state_mut().app.new_experiment_dialog.current_dir =
        path.parent().unwrap().to_path_buf();
    harness.run();
    harness
        .state_mut()
        .app
        .new_experiment_dialog
        .test_next_open_file = Some(path.to_path_buf());
    click_menu_entry(harness, "File Manager", "Import Data");
    harness.get_by_label("Experiment ID:").focus();
    harness.key_press(egui::Key::Home);
    harness.key_press_modifiers(egui::Modifiers::COMMAND, egui::Key::A);
    harness.key_press(egui::Key::Backspace);
    harness
        .get_by_label("Experiment ID:")
        .type_text(experiment_id);
    if let Some(temperature) = isothermal_temperature {
        let input = harness.get_by_label("Isothermal Temperature (K):");
        input.focus();
        harness.key_press_modifiers(egui::Modifiers::COMMAND, egui::Key::A);
        harness.key_press(egui::Key::Backspace);
        input.type_text(&temperature.to_string());
    }
    if let Some(beta) = heating_rate {
        let input = harness.get_by_label("Heating Rate (K/min):");
        input.focus();
        harness.key_press_modifiers(egui::Modifiers::COMMAND, egui::Key::A);
        harness.key_press(egui::Key::Backspace);
        input.type_text(&beta.to_string());
    }
    harness.run();
    harness.get_by_label("Import experiment").click();
    harness.run();
    assert!(
        harness
            .state()
            .app
            .new_experiment_dialog
            .last_error
            .is_none(),
        "GUI import failed: {:?}",
        harness.state().app.new_experiment_dialog.last_error
    );
}

pub fn read_column(model: &PlotModel, id: &str, column: &str) -> Vec<f64> {
    let experiment = model.series.get_experiment_by_id(id).unwrap();
    let frame = experiment.dataset.frame.clone().collect().unwrap();
    frame
        .column(column)
        .unwrap()
        .f64()
        .unwrap()
        .into_no_null_iter()
        .collect()
}

pub fn rmse(actual: &[f64], expected: &[f64]) -> f64 {
    assert_eq!(actual.len(), expected.len(), "RMSE vector length mismatch");
    assert!(!actual.is_empty(), "RMSE requires data");
    (actual
        .iter()
        .zip(expected)
        .map(|(left, right)| (left - right).powi(2))
        .sum::<f64>()
        / actual.len() as f64)
        .sqrt()
}

pub fn ideal_conversion(
    time_seconds: &[f64],
    clean_mass: &[f64],
    baseline_end_seconds: f64,
) -> Vec<f64> {
    assert_eq!(time_seconds.len(), clean_mass.len());
    let baseline = time_seconds
        .iter()
        .zip(clean_mass)
        .filter_map(|(time, mass)| (*time <= baseline_end_seconds).then_some(*mass))
        .collect::<Vec<_>>();
    assert!(!baseline.is_empty(), "conversion baseline has no samples");
    let initial_mass = baseline.iter().sum::<f64>() / baseline.len() as f64;
    clean_mass
        .iter()
        .map(|mass| 1.0 - mass / initial_mass)
        .collect()
}

pub fn quality_from_slices(
    raw: &[f64],
    filtered: &[f64],
    raw_name: &str,
    filtered_name: &str,
) -> FilterQualityReport {
    filter_quality_from_arrays(
        &Array1::from_vec(raw.to_vec()),
        &Array1::from_vec(filtered.to_vec()),
        raw_name,
        filtered_name,
    )
    .expect("compute objective filter quality")
}

pub fn print_filter_report(
    seed: u64,
    experiment_id: &str,
    stage: &str,
    window_points: usize,
    polynomial_order: usize,
    report: &FilterQualityReport,
    truth_rmse: Option<(f64, f64)>,
) {
    if let Some((raw, filtered)) = truth_rmse {
        let improvement = 100.0 * (1.0 - filtered / raw);
        println!(
            "[STORY][filter] seed={seed} experiment={experiment_id} stage={stage} method=Savitzky-Golay window={window_points}points poly={polynomial_order} | {report} | rmse_to_truth raw={raw:.6} filtered={filtered:.6} reduction={improvement:.2}%"
        );
    } else {
        println!(
            "[STORY][filter] seed={seed} experiment={experiment_id} stage={stage} method=Savitzky-Golay window={window_points}points poly={polynomial_order} | {report}"
        );
    }
}

pub fn make_clean_export_series(model: &PlotModel, ids: &[String]) -> TGASeries {
    let columns = ["t", "T", "sg_mass", "eta", "deta_dt"];
    let mut clean = TGASeries::new();
    for id in ids {
        let source = model.series.get_experiment_by_id(id).unwrap();
        let mut experiment = source.clone();
        experiment.dataset.frame = source
            .dataset
            .frame
            .clone()
            .select(columns.iter().map(|name| col(*name)).collect::<Vec<_>>());
        experiment
            .dataset
            .schema
            .columns
            .retain(|name, _| columns.contains(&name.as_str()));
        experiment.dataset.oneframeplot = None;
        experiment.dataset.history_of_operations =
            crate::Kinetics::experimental_kinetics::one_experiment_dataset::History::new();
        experiment.dataset.undo_stack.clear();
        experiment.dataset.initialize_column_provenance();
        clean.push(experiment).unwrap();
    }
    clean
}
