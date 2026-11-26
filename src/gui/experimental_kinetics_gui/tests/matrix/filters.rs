use super::super::support::{read_column, rmse};
use super::generation::{NoiseProfile, ThermalScenario, noisy_and_truth};
use crate::Kinetics::experimental_kinetics::math::exp_kinetics_smooth_filter::HampelStrategy;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::gui::experimental_kinetics_gui::model::{Colours, PlotModel};
use RustedSciThe::numerical::data_processing::lowess_wrapper::{Backend, LowessConfig};
use RustedSciThe::numerical::data_processing::splines::SplineKind;
use tabled::{Table, Tabled, settings::Style};

const SEED_BASE: u64 = 82_000;
const EXPERIMENT_ID: &str = "filter_matrix";

#[derive(Clone, Copy)]
pub(super) enum FilterKind {
    SavitzkyGolay,
    HampelMad,
    RollingMean,
    Lowess,
    LinearResample,
    LsqSpline,
}

impl FilterKind {
    pub(super) const ALL: [Self; 6] = [
        Self::SavitzkyGolay,
        Self::HampelMad,
        Self::RollingMean,
        Self::Lowess,
        Self::LinearResample,
        Self::LsqSpline,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::SavitzkyGolay => "SG 31/3",
            Self::HampelMad => "Hampel-MAD 11/3",
            Self::RollingMean => "Rolling mean 11",
            Self::Lowess => "LOWESS frac=0.12",
            Self::LinearResample => "Linear spline resample",
            Self::LsqSpline => "LSQ spline p=3,k=12",
        }
    }

    pub(super) fn output_column(self) -> &'static str {
        match self {
            Self::SavitzkyGolay => "sg_mass",
            Self::HampelMad => "hampel_mass",
            Self::RollingMean => "rolling_mass",
            Self::Lowess => "lowess_mass",
            Self::LinearResample => "spline_mass",
            Self::LsqSpline => "lsq_mass",
        }
    }
}

#[derive(Tabled)]
struct FilterRow {
    #[tabled(rename = "Seed")]
    seed: u64,
    #[tabled(rename = "Noise profile")]
    noise: String,
    #[tabled(rename = "Filter / parameters")]
    filter: String,
    #[tabled(rename = "r1 ratio")]
    rough1: String,
    #[tabled(rename = "r2 ratio")]
    rough2: String,
    #[tabled(rename = "RMSE input (mg)")]
    rmse_input: String,
    #[tabled(rename = "NRMSE input")]
    nrmse_input: String,
    #[tabled(rename = "Corr input")]
    correlation: String,
    #[tabled(rename = "Truth RMSE raw -> filtered (mg)")]
    truth_rmse: String,
    #[tabled(rename = "Truth error reduction")]
    truth_reduction: String,
}

pub(super) fn model_for_run(
    run: &crate::Kinetics::experimental_kinetics::testing_mod::SimulatedTGA,
) -> PlotModel {
    let mut model = PlotModel::new();
    let mut meta = run.meta.clone();
    meta.id = EXPERIMENT_ID.to_string();
    model
        .create_from_synthetic_data(&run.virtual_tga, meta)
        .expect("create test experiment");
    model
        .series
        .bind_time(EXPERIMENT_ID, "time", Unit::Second)
        .unwrap();
    model
        .series
        .bind_temperature(EXPERIMENT_ID, "temperature", Unit::Kelvin)
        .unwrap();
    model
        .series
        .bind_mass(EXPERIMENT_ID, "mass", Unit::Milligram)
        .unwrap();
    model.set_x(EXPERIMENT_ID, "time").unwrap();
    model.set_y(EXPERIMENT_ID, "mass").unwrap();
    model
        .create_points_for_curve_with_builder(EXPERIMENT_ID, Colours::Blue)
        .unwrap();
    model.select_curve(0);
    model
}

pub(super) fn apply_filter(model: &mut PlotModel, filter: FilterKind, n_points: usize) -> String {
    let out = filter.output_column();
    match filter {
        FilterKind::SavitzkyGolay => model
            .sg_filter_column_for_selected_as("mass", 31, 3, 0, 1.0, Some(out))
            .expect("GUI model handler: Savitzky-Golay"),
        FilterKind::HampelMad => model
            .hampel_filter_for_selected_as(
                "mass",
                11,
                3.0,
                HampelStrategy::ReplaceWithMedian,
                Some(out),
            )
            .expect("GUI model handler: Hampel/MAD"),
        FilterKind::RollingMean => model
            .rolling_mean_for_selected_as("mass", 11, Some(out))
            .expect("GUI model handler: rolling mean"),
        FilterKind::Lowess => model
            .lowess_smooth_columns_for_selected_as(
                "time",
                &["mass"],
                &[Some(out)],
                LowessConfig {
                    fraction: 0.12,
                    iterations: 2,
                    delta: 0.0,
                    parallel: false,
                    backend: Backend::CPUBack,
                },
            )
            .expect("GUI model handler: LOWESS"),
        FilterKind::LinearResample => model
            .splines_for_selected_as("t_spline", n_points, SplineKind::Linear, Some(out))
            .expect("GUI model handler: linear spline resampling"),
        FilterKind::LsqSpline => model
            .lsq_spline_resample_columns_as(
                "time",
                "t_lsq",
                &["mass"],
                &[Some(out)],
                n_points,
                3,
                12,
            )
            .expect("GUI model handler: LSQ spline"),
    }
    out.to_string()
}

pub(super) fn truth_at(time: &[f64], truth_time: &[f64], truth_mass: &[f64]) -> Vec<f64> {
    assert_eq!(truth_time.len(), truth_mass.len());
    time.iter()
        .map(|&x| {
            let upper = truth_time.partition_point(|&value| value < x);
            if upper == 0 {
                return truth_mass[0];
            }
            if upper == truth_time.len() {
                return truth_mass[truth_mass.len() - 1];
            }
            let lo = upper - 1;
            let fraction = (x - truth_time[lo]) / (truth_time[upper] - truth_time[lo]);
            truth_mass[lo] + fraction * (truth_mass[upper] - truth_mass[lo])
        })
        .collect()
}

#[test]
fn gui_filter_stage_compares_supported_smoothers_across_noise_profiles() {
    for scenario in [ThermalScenario::Isothermal, ThermalScenario::NonIsothermal] {
        let mut rows = Vec::new();
        for (profile_index, noise) in NoiseProfile::ALL.into_iter().enumerate() {
            let seed = SEED_BASE + profile_index as u64 * 31 + scenario as u64;
            let (runs, truth) = noisy_and_truth(scenario, noise, seed);
            let noisy_run = &runs[0];
            let clean_run = &truth[0];
            let baseline_model = model_for_run(noisy_run);
            let mut app = super::super::support::empty_app();
            app.model = baseline_model.clone();

            super::super::support::with_app(app, |harness| {
                for filter in FilterKind::ALL {
                    harness.state_mut().app.model = baseline_model.clone();
                    harness.run();
                    let n_points = noisy_run.virtual_tga.time.len();
                    let out = apply_filter(&mut harness.state_mut().app.model, filter, n_points);
                    harness.run();

                    let model = &harness.state().app.model;
                    let dataset = &model
                        .series
                        .get_experiment_by_id(EXPERIMENT_ID)
                        .unwrap()
                        .dataset;
                    let frame = dataset.frame.clone().collect().unwrap();
                    assert!(
                        frame.column(&out).is_ok(),
                        "{} did not create {out}",
                        filter.label()
                    );
                    assert!(
                        dataset.column_provenance_text(&out).is_some(),
                        "{} did not record provenance",
                        filter.label()
                    );

                    let paired_samples = dataset
                        .preview_rows(&["mass", &out], frame.height())
                        .expect("read aligned nullable filter samples");
                    let (paired_raw, paired_filtered): (Vec<_>, Vec<_>) =
                        paired_samples.columns["mass"]
                            .iter()
                            .zip(&paired_samples.columns[&out])
                            .filter_map(|(raw, filtered)| match (raw, filtered) {
                                (Some(raw), Some(filtered)) => Some((*raw, *filtered)),
                                _ => None,
                            })
                            .unzip();
                    let paired_raw = ndarray::Array1::from_vec(paired_raw);
                    let paired_filtered = ndarray::Array1::from_vec(paired_filtered);
                    let report = crate::Kinetics::experimental_kinetics::math::ndarray_statistics::filter_quality_from_arrays(
                        &paired_raw,
                        &paired_filtered,
                        "mass",
                        &out,
                    )
                    .expect("filter quality report on shared valid rows");
                    let raw_time = read_column(model, EXPERIMENT_ID, "time");
                    let raw_mass = read_column(model, EXPERIMENT_ID, "mass");
                    let filtered_time_col = match filter {
                        FilterKind::LinearResample => "t_spline",
                        FilterKind::LsqSpline => "t_lsq",
                        _ => "time",
                    };
                    let filtered_samples = dataset
                        .preview_rows(&[filtered_time_col, &out], frame.height())
                        .expect("read filtered values on their shared valid rows");
                    let filtered_time_values = &filtered_samples.columns[filtered_time_col];
                    let (filtered_time, filtered_mass): (Vec<_>, Vec<_>) = filtered_time_values
                        .iter()
                        .zip(&filtered_samples.columns[&out])
                        .filter_map(|(time, mass)| match (time, mass) {
                            (Some(time), Some(mass)) => Some((*time, *mass)),
                            _ => None,
                        })
                        .unzip();
                    assert_eq!(filtered_time.len(), filtered_mass.len());
                    let raw_truth = truth_at(
                        &raw_time,
                        &clean_run.virtual_tga.time,
                        &clean_run.virtual_tga.mass,
                    );
                    let filtered_truth = truth_at(
                        &filtered_time,
                        &clean_run.virtual_tga.time,
                        &clean_run.virtual_tga.mass,
                    );
                    let raw_rmse = rmse(&raw_mass, &raw_truth);
                    let filtered_rmse = rmse(&filtered_mass, &filtered_truth);
                    let reduction = 100.0 * (1.0 - filtered_rmse / raw_rmse);
                    assert!(
                        report.roughness_ratio_1st.is_finite()
                            && report.roughness_ratio_2nd.is_finite()
                            && report.rmse.is_finite()
                            && report.normalized_rmse.is_finite()
                            && report.correlation.is_finite(),
                        "non-finite report for {} / {}: {report}",
                        noise.label(),
                        filter.label()
                    );
                    rows.push(FilterRow {
                        seed,
                        noise: noise.label().into(),
                        filter: filter.label().into(),
                        rough1: format!("{:.3}", report.roughness_ratio_1st),
                        rough2: format!("{:.3}", report.roughness_ratio_2nd),
                        rmse_input: format!("{:.5}", report.rmse),
                        nrmse_input: format!("{:.4}", report.normalized_rmse),
                        correlation: format!("{:.5}", report.correlation),
                        truth_rmse: format!("{raw_rmse:.4} -> {filtered_rmse:.4}"),
                        truth_reduction: format!("{reduction:.1}%"),
                    });
                }
            });
        }
        println!(
            "\n=== STORY MATRIX / 2. GUI FILTERS / {} ===\n{}",
            scenario.label(),
            Table::new(rows).with(Style::rounded())
        );
    }
}
