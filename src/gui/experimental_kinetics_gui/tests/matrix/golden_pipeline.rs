use super::super::support::{empty_app, quality_from_slices, read_column, rmse, with_app};
use super::filters::{model_for_run, truth_at};
use super::generation::{NoiseProfile, ThermalScenario, noisy_and_truth};
use crate::Kinetics::experimental_kinetics::exp_engine_api::{
    GoldenPipelineConfig, SavGolConfig, SplineConfig,
};
use egui_kittest::kittest::Queryable;
use tabled::{Table, Tabled, settings::Style};

const SEED_BASE: u64 = 96_000;
const OUTPUT_ID: &str = "proceeded_filter_matrix";

#[derive(Tabled)]
struct GoldenPipelineRow {
    #[tabled(rename = "Seed")]
    seed: u64,
    #[tabled(rename = "Thermal case")]
    scenario: String,
    #[tabled(rename = "Noise profile")]
    noise: String,
    #[tabled(rename = "Golden SG window / poly")]
    sg: String,
    #[tabled(rename = "Samples")]
    samples: String,
    #[tabled(rename = "rough1")]
    rough1: String,
    #[tabled(rename = "rough2")]
    rough2: String,
    #[tabled(rename = "RMSE vs raw input (mg)")]
    input_rmse: String,
    #[tabled(rename = "NRMSE vs SG input")]
    input_nrmse: String,
    #[tabled(rename = "Correlation vs SG input")]
    input_correlation: String,
    #[tabled(rename = "Truth RMSE input -> output, shared grid (mg)")]
    truth_rmse: String,
    #[tabled(rename = "Truth error reduction")]
    truth_reduction: String,
    #[tabled(rename = "GUI outcome")]
    outcome: String,
}

fn golden_config(window: usize, input_rows: usize) -> GoldenPipelineConfig {
    GoldenPipelineConfig {
        k: 1.0,
        b: 0.0,
        time_cut_before: None,
        time_cut_after: None,
        sav_gol_config: SavGolConfig {
            window_size: window,
            poly_degree: 3,
            deriv: 0,
            delta: 1.0,
        },
        averaging_time: 10.0,
        spline_config: SplineConfig {
            n_points: input_rows * 2 / 3,
            n_internal_points: 12,
            degree: 3,
        },
        save_to_new_experiment: true,
        del_old_experiment: true,
    }
}

fn run_golden_pipeline_from_gui(
    app: crate::gui::experimental_kinetics_gui::experimental_kinetics_gui_main::PlotApp,
    config: GoldenPipelineConfig,
    scenario: ThermalScenario,
    noise: NoiseProfile,
    seed: u64,
    noisy_time: &[f64],
    noisy_mass: &[f64],
    clean_time: &[f64],
    clean_mass: &[f64],
    window: usize,
) -> Result<GoldenPipelineRow, String> {
    let mut outcome = Err("Golden Pipeline Calculate was not applied".to_string());
    with_app(app, |harness| {
        harness.set_size(eframe::egui::vec2(1600.0, 1200.0));
        harness.run();
        harness.get_by_label("golden pipeline").click();
        harness.run();
        assert!(
            harness
                .state()
                .app
                .quick_actions_state
                .golden_pipeline_dialog
                .open,
            "front-panel action did not open the Golden Pipeline dialog"
        );
        assert!(
            harness.query_by_label("Calculate").is_some(),
            "Golden Pipeline dialog did not expose its Calculate action"
        );
        harness
            .state_mut()
            .app
            .quick_actions_state
            .golden_pipeline_dialog
            .config = config;
        harness.run();
        harness.get_by_label("Calculate").click();
        harness.run();

        let model = &harness.state().app.model;
        let dialog = &harness
            .state()
            .app
            .quick_actions_state
            .golden_pipeline_dialog;
        if let Some(error) = &dialog.last_error {
            outcome = Err(error.clone());
            return;
        }
        let Some(experiment) = model.series.get_experiment_by_id(OUTPUT_ID).ok() else {
            outcome = Err(format!(
                "Calculate returned without an error but did not create '{OUTPUT_ID}'"
            ));
            return;
        };
        let Some(time_col) = experiment.dataset.schema.time.as_deref() else {
            outcome = Err("Golden Pipeline output has no bound time column".into());
            return;
        };
        let Some(mass_col) = experiment.dataset.schema.mass.as_deref() else {
            outcome = Err("Golden Pipeline output has no bound mass column".into());
            return;
        };

        let output_time_hours = read_column(model, OUTPUT_ID, time_col);
        let output_mass = read_column(model, OUTPUT_ID, mass_col);
        if output_time_hours.len() != output_mass.len() || output_mass.is_empty() {
            outcome = Err(format!(
                "Output time/mass lengths are invalid: {} vs {}",
                output_time_hours.len(),
                output_mass.len()
            ));
            return;
        }
        let output_time_seconds = output_time_hours
            .iter()
            .map(|time| time * 3600.0)
            .collect::<Vec<_>>();
        let noisy_at_output = truth_at(&output_time_seconds, noisy_time, noisy_mass);
        let clean_at_output = truth_at(&output_time_seconds, clean_time, clean_mass);
        let raw_truth_rmse = rmse(&noisy_at_output, &clean_at_output);
        let output_truth_rmse = rmse(&output_mass, &clean_at_output);
        let reduction = 100.0 * (1.0 - output_truth_rmse / raw_truth_rmse);
        let report = quality_from_slices(
            &noisy_at_output,
            &output_mass,
            "input_mass_on_output_grid",
            "golden_pipeline_mass",
        );
        if ![
            report.roughness_ratio_1st,
            report.roughness_ratio_2nd,
            report.rmse,
            report.normalized_rmse,
            report.correlation,
            output_truth_rmse,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            outcome = Err("Golden Pipeline produced non-finite quality metrics".into());
            return;
        }
        if !dialog.open {
            outcome = Err("Golden Pipeline closed before showing its quality report".into());
            return;
        }
        if harness
            .query_by_label("Savitzky-Golay mass smoothing quality")
            .is_none()
        {
            outcome = Err("Golden Pipeline quality summary is not visible in the dialog".into());
            return;
        }
        for metric in [
            "First-difference roughness ratio",
            "Second-difference roughness ratio",
            "RMSE vs input (mass units)",
            "Normalized RMSE",
            "Correlation",
        ] {
            if harness.query_by_label(metric).is_none() {
                outcome = Err(format!("Golden Pipeline report is missing '{metric}'"));
                return;
            }
        }
        let Some(ui_report) = dialog.quality_report.as_ref() else {
            outcome = Err("Golden Pipeline dialog did not retain its quality report".into());
            return;
        };
        for value in [
            ui_report.roughness_ratio_1st,
            ui_report.roughness_ratio_2nd,
            ui_report.rmse,
            ui_report.normalized_rmse,
            ui_report.correlation,
        ] {
            let formatted_value = format!("{value:.6}");
            if harness.query_by_label(&formatted_value).is_none() {
                outcome = Err(format!(
                    "Golden Pipeline dialog does not show metric {value:.6}"
                ));
                return;
            }
        }

        outcome = Ok(GoldenPipelineRow {
            seed,
            scenario: scenario.label().into(),
            noise: noise.label().into(),
            sg: format!("{window}/3"),
            samples: output_mass.len().to_string(),
            rough1: format!("{:.3}", report.roughness_ratio_1st),
            rough2: format!("{:.3}", report.roughness_ratio_2nd),
            input_rmse: format!("{:.5}", report.rmse),
            input_nrmse: format!("{:.5}", report.normalized_rmse),
            input_correlation: format!("{:.5}", report.correlation),
            truth_rmse: format!("{raw_truth_rmse:.4} -> {output_truth_rmse:.4}"),
            truth_reduction: format!("{reduction:.1}%"),
            outcome: format!("success; output={OUTPUT_ID}"),
        });
    });
    outcome
}

#[test]
fn golden_pipeline_dialog_runs_quantitative_filtering_matrix() {
    let windows = [11, 31, 51];
    let mut rows = Vec::new();
    let mut failures = Vec::new();

    for scenario in [ThermalScenario::Isothermal, ThermalScenario::NonIsothermal] {
        for (noise_index, noise) in NoiseProfile::ALL.into_iter().enumerate() {
            let seed = SEED_BASE + (scenario as u64) * 100 + noise_index as u64;
            let (noisy_runs, clean_runs) = noisy_and_truth(scenario, noise, seed);
            let noisy = &noisy_runs[0];
            let clean = &clean_runs[0];
            for window in windows {
                let baseline_model = model_for_run(noisy);
                let mut app = empty_app();
                app.model = baseline_model;
                app.model.settings.set_calibration_coeffs(1.0, 0.0);
                let result = run_golden_pipeline_from_gui(
                    app,
                    golden_config(window, noisy.virtual_tga.time.len()),
                    scenario,
                    noise,
                    seed,
                    &noisy.virtual_tga.time,
                    &noisy.virtual_tga.mass,
                    &clean.virtual_tga.time,
                    &clean.virtual_tga.mass,
                    window,
                );
                match result {
                    Ok(row) => rows.push(row),
                    Err(error) => {
                        failures.push(format!(
                            "{} / {} / window={window}: {error}",
                            scenario.label(),
                            noise.label()
                        ));
                        rows.push(GoldenPipelineRow {
                            seed,
                            scenario: scenario.label().into(),
                            noise: noise.label().into(),
                            sg: format!("{window}/3"),
                            samples: "-".into(),
                            rough1: "-".into(),
                            rough2: "-".into(),
                            input_rmse: "-".into(),
                            input_nrmse: "-".into(),
                            input_correlation: "-".into(),
                            truth_rmse: "-".into(),
                            truth_reduction: "-".into(),
                            outcome: format!("FAILED: {error}"),
                        });
                    }
                }
            }
        }
    }

    println!(
        "\n=== STORY MATRIX / GOLDEN PIPELINE GUI (SG window sweep; paired clean truth) ===\n{}",
        Table::new(rows).with(Style::rounded())
    );
    assert!(
        failures.is_empty(),
        "Golden Pipeline GUI failures:\n{}",
        failures.join("\n")
    );
}
