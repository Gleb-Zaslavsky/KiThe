use super::super::support::{empty_app, read_column, with_app};
use super::filters::{FilterKind, apply_filter};
use super::generation::{NoiseProfile, ThermalScenario, config, noisy_and_truth};
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::Kinetics::experimental_kinetics::kinetic_methods::isoconversion::IsoconversionalMethod;
use crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::linear_regression;
use crate::Kinetics::experimental_kinetics::math::fitting::FittingModelName;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use crate::Kinetics::experimental_kinetics::testing_mod::{
    ExperimentMode, SimulatedTGA, VirtualTGA,
};
use crate::gui::experimental_kinetics_gui::model::{Colours, PlotModel};
use egui_kittest::kittest::Queryable;
use ndarray::Array1;
use tabled::{Table, Tabled, settings::Style};

const R_GAS: f64 = 8.314;
const E_TRUE: f64 = 85_000.0;
const A_TRUE: f64 = 2.0e5;

#[derive(Tabled)]
struct KineticRow {
    #[tabled(rename = "Thermal case")]
    scenario: String,
    #[tabled(rename = "Noise / seed")]
    noise_seed: String,
    #[tabled(rename = "Mass filter")]
    filter: String,
    #[tabled(rename = "Kinetic method")]
    method: String,
    #[tabled(rename = "Fit / method R2")]
    r2: String,
    #[tabled(rename = "Ea (kJ/mol), true=85")]
    ea: String,
    #[tabled(rename = "Ea error")]
    ea_error: String,
    #[tabled(rename = "A (1/s), true=2e5")]
    pre_exponential: String,
    #[tabled(rename = "A error")]
    a_error: String,
    #[tabled(rename = "Outcome")]
    outcome: String,
}

fn model_for_runs(runs: &[SimulatedTGA]) -> (PlotModel, Vec<String>) {
    let mut model = PlotModel::new();
    let mut ids = Vec::with_capacity(runs.len());
    for (index, run) in runs.iter().enumerate() {
        let id = format!("matrix_run_{index}");
        let meta = ExperimentMeta {
            id: id.clone(),
            heating_rate: run.meta.heating_rate,
            isothermal_temperature: run.meta.isothermal_temperature,
            thermal_mode: run.meta.thermal_mode,
            comment: run.meta.comment.clone(),
        };
        model
            .create_from_synthetic_data(&run.virtual_tga, meta)
            .expect("create matrix experiment");
        model.series.bind_time(&id, "time", Unit::Second).unwrap();
        model
            .series
            .bind_temperature(&id, "temperature", Unit::Kelvin)
            .unwrap();
        model
            .series
            .bind_mass(&id, "mass", Unit::Milligram)
            .unwrap();
        model.set_x(&id, "time").unwrap();
        model.set_y(&id, "mass").unwrap();
        model
            .create_points_for_curve_with_builder(&id, Colours::Blue)
            .unwrap();
        ids.push(id);
    }
    (model, ids)
}

fn select_experiment_curve(model: &mut PlotModel, id: &str) {
    let index = model
        .plots
        .iter()
        .position(|curve| curve.experiment_id == id)
        .expect("each matrix experiment has a visible GUI curve");
    model.select_curve(index);
}

fn run_fit_gui(
    base: &PlotModel,
    id: &str,
    y_col: &str,
    output: &str,
    model_name: FittingModelName,
) -> Result<(f64, f64), String> {
    let mut app = empty_app();
    app.model = base.clone();
    let mut result = Err("Fit window did not run".to_string());
    with_app(app, |harness| {
        super::super::support::click_menu_entry(harness, "Kinetic Methods", "Fit Model");
        {
            let state = &mut harness.state_mut().app.kinetic_methods_state;
            state.fit_model_experiment_id = id.to_string();
            state.fit_model_x_col = "time".to_string();
            state.fit_model_y_col = y_col.to_string();
            state.fit_model_output_col = output.to_string();
            state.fit_model_model = model_name;
            state.fit_model_initial_guess.clear();
        }
        harness.run();
        harness.get_by_label("Fit").click();
        harness.run();
        result = match harness
            .state()
            .app
            .kinetic_methods_state
            .fit_model_result
            .as_ref()
        {
            Some(fit) => fit
                .coefficients
                .iter()
                .find_map(|(name, value)| (name == "k").then_some((*value, fit.r2)))
                .ok_or_else(|| "fit result has no k coefficient".to_string()),
            None => Err(harness
                .state()
                .app
                .kinetic_methods_state
                .fit_model_status
                .clone()),
        };
    });
    result
}

fn run_isoconversional_gui(
    base: &PlotModel,
    ids: &[String],
    method: IsoconversionalMethod,
    output_id: &str,
) -> Result<(f64, f64, f64, usize), String> {
    let mut app = empty_app();
    app.model = base.clone();
    let mut result = Err("Isoconversional window did not run".to_string());
    with_app(app, |harness| {
        {
            let state = &mut harness.state_mut().app.kinetic_methods_state;
            state.isoconversional_method = method;
            state.isoconversional_output_id = output_id.to_string();
            state.isoconversional_use_selected_experiments = true;
            state.isoconversional_selected_experiments = ids.to_vec();
        }
        harness.run();
        super::super::support::click_menu_entry(harness, "Kinetic Methods", "Isoconversional");
        harness.get_by_label("Calculate").click();
        harness.run();
        let state = &harness.state().app.kinetic_methods_state;
        if !state.isoconversional_status.starts_with("Done.") {
            result = Err(state.isoconversional_status.clone());
            return;
        }
        let ea = read_column(&harness.state().app.model, output_id, "Ea");
        let alpha = read_column(&harness.state().app.model, output_id, "eta");
        let r2 = read_column(&harness.state().app.model, output_id, "R2");
        let values = ea
            .into_iter()
            .zip(alpha)
            .zip(r2)
            .filter_map(|((value, conversion), r2)| {
                (value.is_finite() && r2.is_finite() && (0.1..=0.9).contains(&conversion))
                    .then_some((value, r2))
            })
            .collect::<Vec<_>>();
        if values.is_empty() {
            result = Err("method returned no finite Ea in conversion range 0.1..0.9".into());
        } else {
            let mean = values.iter().map(|(value, _)| value).sum::<f64>() / values.len() as f64;
            let variance = values
                .iter()
                .map(|(value, _)| (value - mean).powi(2))
                .sum::<f64>()
                / values.len() as f64;
            let mean_r2 = values.iter().map(|(_, r2)| r2).sum::<f64>() / values.len() as f64;
            result = Ok((mean, variance.sqrt(), mean_r2, values.len()));
        }
    });
    result
}

fn arrhenius_regression(temperatures: &[f64], rates_per_second: &[f64]) -> (f64, f64, f64) {
    let regression = linear_regression(
        &Array1::from_vec(temperatures.iter().map(|temp| 1.0 / temp).collect()),
        &Array1::from_vec(rates_per_second.iter().map(|rate| rate.ln()).collect()),
    )
    .expect("fit Arrhenius parameters from finite positive temperatures and rates");
    (
        -regression.slope * R_GAS,
        regression.intercept.exp(),
        regression.r2,
    )
}

fn row_for_result(
    scenario: ThermalScenario,
    noise: NoiseProfile,
    seed: u64,
    filter: FilterKind,
    method: &str,
    r2: Option<f64>,
    ea_j_mol: Option<f64>,
    ea_sd_j_mol: Option<f64>,
    a: Option<f64>,
    outcome: &str,
) -> KineticRow {
    let ea = ea_j_mol
        .map(|value| format!("{:.3}", value / 1000.0))
        .unwrap_or_else(|| "-".into());
    let ea_error = ea_j_mol
        .map(|value| format!("{:.1}%", 100.0 * (value - E_TRUE).abs() / E_TRUE))
        .unwrap_or_else(|| "-".into());
    let a_error = a
        .map(|value| format!("{:.1}%", 100.0 * ((value / A_TRUE).ln()).exp_m1().abs()))
        .unwrap_or_else(|| "-".into());
    KineticRow {
        scenario: scenario.label().into(),
        noise_seed: format!("{} / {seed}", noise.label()),
        filter: filter.label().into(),
        method: method.into(),
        r2: r2
            .map(|value| format!("{value:.5}"))
            .unwrap_or_else(|| "-".into()),
        ea: match (ea_sd_j_mol, ea.as_str()) {
            (Some(sd), value) if value != "-" => format!("{value} +/- {:.3}", sd / 1000.0),
            _ => ea,
        },
        ea_error,
        pre_exponential: a
            .map(|value| format!("{value:.4e}"))
            .unwrap_or_else(|| "n/a".into()),
        a_error,
        outcome: outcome.into(),
    }
}

#[test]
fn gui_kinetic_stage_prints_isothermal_and_nonisothermal_method_matrix() {
    let mut rows = Vec::new();
    let cases = [
        (NoiseProfile::Gaussian, FilterKind::SavitzkyGolay),
        (NoiseProfile::Impulse, FilterKind::HampelMad),
        (NoiseProfile::Flicker, FilterKind::Lowess),
    ];

    for (case_index, (noise, filter)) in cases.into_iter().enumerate() {
        let isothermal_seed = 93_000 + case_index as u64 * 17;
        let (isothermal_runs, _) =
            noisy_and_truth(ThermalScenario::Isothermal, noise, isothermal_seed);
        let (mut model, iso_ids) = model_for_runs(&isothermal_runs);
        let mut app = empty_app();
        app.model = model.clone();
        with_app(app, |harness| {
            for id in &iso_ids {
                select_experiment_curve(&mut harness.state_mut().app.model, id);
                harness.run();
                let out = apply_filter(
                    &mut harness.state_mut().app.model,
                    filter,
                    isothermal_runs[0].virtual_tga.time.len(),
                );
                harness.run();
                let current = &mut harness.state_mut().app.model;
                current.series.bind_mass(id, &out, Unit::Milligram).unwrap();
                current.from_s_to_h_for_experiment(id).unwrap();
                current
                    .conversion_for_experiment(id.clone(), 0.005, "eta")
                    .unwrap();
                current
                    .series
                    .sg_filter_column_as(id, "eta", 31, 3, 0, 1.0, Some("eta"))
                    .unwrap();
                current.derive_deta_dt(id, Some("deta_dt")).unwrap();
                current.series.drop_nulls(id).unwrap();
                current
                    .series
                    .sg_filter_column_as(id, "deta_dt", 11, 3, 0, 1.0, Some("deta_dt"))
                    .unwrap();
            }
            model = harness.state().app.model.clone();
        });

        let mut temperatures = Vec::new();
        let mut rates = Vec::new();
        let mut fit_r2 = Vec::new();
        let mut fit_failures = Vec::new();
        for (id, run) in iso_ids.iter().zip(&isothermal_runs) {
            let temp = run
                .meta
                .isothermal_temperature
                .expect("isothermal metadata");
            match run_fit_gui(
                &model,
                id,
                filter.output_column(),
                &format!("fit_{id}"),
                FittingModelName::DecExp,
            ) {
                Ok((k_per_hour, r2)) if k_per_hour > 0.0 && k_per_hour.is_finite() => {
                    temperatures.push(temp);
                    rates.push(k_per_hour / 3600.0);
                    fit_r2.push(r2);
                }
                Ok((k, r2)) => fit_failures.push(format!("{id}: invalid k={k}, R2={r2}")),
                Err(error) => fit_failures.push(format!("{id}: {error}")),
            }
        }
        if temperatures.len() >= 2 {
            let (ea, a, r2) = arrhenius_regression(&temperatures, &rates);
            let mean_fit_r2 = fit_r2.iter().sum::<f64>() / fit_r2.len() as f64;
            rows.push(row_for_result(
                ThermalScenario::Isothermal,
                noise,
                isothermal_seed,
                filter,
                "GUI DecExp + Arrhenius",
                Some(r2),
                Some(ea),
                None,
                Some(a),
                &format!(
                    "{} of 3 fits; mean DecExp R2={mean_fit_r2:.5}",
                    fit_r2.len()
                ),
            ));
        } else {
            rows.push(row_for_result(
                ThermalScenario::Isothermal,
                noise,
                isothermal_seed,
                filter,
                "GUI DecExp + Arrhenius",
                None,
                None,
                None,
                None,
                &fit_failures.join("; "),
            ));
        }

        // Friedman Integral is the isoconversional method configured for isothermal data.
        for (offset, method) in [IsoconversionalMethod::FriedmanIntegral]
            .into_iter()
            .enumerate()
        {
            let output = format!("iso_result_{case_index}_{offset}");
            match run_isoconversional_gui(&model, &iso_ids, method, &output) {
                Ok((mean_ea, sd_ea, mean_r2, count)) => rows.push(row_for_result(
                    ThermalScenario::Isothermal,
                    noise,
                    isothermal_seed,
                    filter,
                    method.display_name(),
                    Some(mean_r2),
                    Some(mean_ea),
                    Some(sd_ea),
                    None,
                    &format!("{count} conversion levels"),
                )),
                Err(error) => rows.push(row_for_result(
                    ThermalScenario::Isothermal,
                    noise,
                    isothermal_seed,
                    filter,
                    method.display_name(),
                    None,
                    None,
                    None,
                    None,
                    &error,
                )),
            }
        }

        let noniso_seed = 94_000 + case_index as u64 * 17;
        let mut noniso_config = config(ThermalScenario::NonIsothermal, noise, noniso_seed);
        noniso_config.n_points = 3_000;
        noniso_config.experiment_mode =
            crate::Kinetics::experimental_kinetics::testing_mod::ExperimentMode::NonIsothermal {
                t0: 350.0,
                heating_rates: vec![2.0, 5.0, 10.0],
            };
        let runs = VirtualTGA::generate_series(&noniso_config).unwrap();
        let (mut noniso_model, noniso_ids) = model_for_runs(&runs);
        let mut app = empty_app();
        app.model = noniso_model.clone();
        with_app(app, |harness| {
            for id in &noniso_ids {
                select_experiment_curve(&mut harness.state_mut().app.model, id);
                harness.run();
                let out = apply_filter(
                    &mut harness.state_mut().app.model,
                    filter,
                    runs[0].virtual_tga.time.len(),
                );
                harness.run();
                let current = &mut harness.state_mut().app.model;
                current.series.bind_mass(id, &out, Unit::Milligram).unwrap();
                current.from_s_to_h_for_experiment(id).unwrap();
                current
                    .conversion_for_experiment(id.clone(), 0.005, "eta")
                    .unwrap();
                current
                    .series
                    .sg_filter_column_as(id, "eta", 31, 3, 0, 1.0, Some("eta"))
                    .unwrap();
                current.derive_deta_dt(id, Some("deta_dt")).unwrap();
                current.series.drop_nulls(id).unwrap();
                current
                    .series
                    .sg_filter_column_as(id, "deta_dt", 11, 3, 0, 1.0, Some("deta_dt"))
                    .unwrap();
            }
            noniso_model = harness.state().app.model.clone();
        });

        let mut branches = vec![("noniso raw-rate validation", None, None, None, None)];
        let mut temperatures = Vec::new();
        let mut rates = Vec::new();
        for id in &noniso_ids {
            let temp = read_column(&noniso_model, id, "temperature");
            let eta = read_column(&noniso_model, id, "eta");
            let rate_h = read_column(&noniso_model, id, "deta_dt");
            for ((temp, eta), rate_h) in temp.into_iter().zip(eta).zip(rate_h) {
                if (0.1..=0.9).contains(&eta) && rate_h > 0.0 {
                    let rate_s = rate_h / 3600.0 / (1.0 - eta);
                    if rate_s.is_finite() && rate_s > 0.0 {
                        temperatures.push(temp);
                        rates.push(rate_s);
                    }
                }
            }
        }
        if temperatures.len() >= 3 {
            let (ea, a, r2) = arrhenius_regression(&temperatures, &rates);
            branches[0] = (
                "Arrhenius check from GUI-processed rate",
                Some(r2),
                Some(ea),
                None,
                Some(a),
            );
        }
        rows.push(row_for_result(
            ThermalScenario::NonIsothermal,
            noise,
            noniso_seed,
            filter,
            branches[0].0,
            branches[0].1,
            branches[0].2,
            branches[0].3,
            branches[0].4,
            "independent validation; not isoconversional output",
        ));

        for (method_index, method) in [
            IsoconversionalMethod::OFW,
            IsoconversionalMethod::KAS,
            IsoconversionalMethod::Starink,
            IsoconversionalMethod::FriedmanDifferential,
            IsoconversionalMethod::Vyazovkin,
        ]
        .into_iter()
        .enumerate()
        {
            let output = format!("noniso_result_{case_index}_{method_index}");
            match run_isoconversional_gui(&noniso_model, &noniso_ids, method, &output) {
                Ok((mean_ea, sd_ea, mean_r2, count)) => {
                    let relative_error = (mean_ea - E_TRUE).abs() / E_TRUE;
                    assert!(
                        relative_error <= 0.12,
                        "{} failed the non-isothermal recovery gate for {} / {}: mean Ea={mean_ea:.1} J/mol, truth={E_TRUE:.1}, relative error={:.2}%",
                        method.display_name(),
                        noise.label(),
                        filter.label(),
                        100.0 * relative_error
                    );
                    rows.push(row_for_result(
                        ThermalScenario::NonIsothermal,
                        noise,
                        noniso_seed,
                        filter,
                        method.display_name(),
                        Some(mean_r2),
                        Some(mean_ea),
                        Some(sd_ea),
                        None,
                        &format!("{count} conversion levels; Ea recovery gate passed"),
                    ));
                }
                Err(error) => {
                    assert!(
                        error.contains("conversion must be strictly increasing"),
                        "{} had an unexpected error for {} / {}: {error}",
                        method.display_name(),
                        noise.label(),
                        filter.label()
                    );
                    rows.push(row_for_result(
                        ThermalScenario::NonIsothermal,
                        noise,
                        noniso_seed,
                        filter,
                        method.display_name(),
                        None,
                        None,
                        None,
                        None,
                        &format!("expected precondition rejection: {error}"),
                    ));
                }
            }
        }
    }

    println!(
        "\n=== STORY MATRIX / 3. GUI KINETIC METHODS (truth E={E_TRUE:.0} J/mol, A={A_TRUE:.3e} 1/s) ===\n{}",
        Table::new(rows).with(Style::rounded())
    );
}

#[test]
fn clean_cold_start_nonisothermal_methods_recover_activation_energy_gate() {
    let seed = 95_001;
    let mut clean_config = config(ThermalScenario::NonIsothermal, NoiseProfile::Gaussian, seed);
    clean_config.n_points = 3_000;
    clean_config.mass_noise = None;
    clean_config.experiment_mode = ExperimentMode::NonIsothermal {
        t0: 350.0,
        heating_rates: vec![2.0, 5.0, 10.0],
    };
    let runs = VirtualTGA::generate_series(&clean_config).expect("generate clean cold-start runs");
    let (mut model, ids) = model_for_runs(&runs);

    // Exercise the same domain preparation used by the GUI. Keep the signal
    // noise-free so this gate measures kinetic recovery, not filter behavior.
    for id in &ids {
        model
            .series
            .sg_filter_column_as(id, "mass", 31, 3, 0, 1.0, Some("sg_mass"))
            .expect("smooth mass on the same path as the GUI story");
        model
            .series
            .bind_mass(id, "sg_mass", Unit::Milligram)
            .expect("bind smoothed mass for conversion");
        model
            .from_s_to_h_for_experiment(id)
            .expect("convert time to hours");
        model
            .conversion_for_experiment(id.clone(), 0.005, "eta")
            .expect("derive conversion from the active mass column");
        let eta = read_column(&model, id, "eta");
        if let Some(index) = eta.iter().position(|value| *value >= 0.95) {
            let time = read_column(&model, id, "time");
            model
                .series
                .try_transform_by_id(id, |experiment| experiment.cut_after_time(time[index]))
                .expect("trim the saturated conversion tail before kinetic analysis");
        }
        model
            .series
            .sg_filter_column_as(id, "eta", 31, 3, 0, 1.0, Some("eta"))
            .expect("smooth conversion on the same path as the GUI story");
        model
            .derive_deta_dt(id, Some("deta_dt"))
            .expect("derive conversion rate");
        model
            .series
            .drop_nulls(id)
            .expect("remove derivative edge rows");
        model
            .series
            .sg_filter_column_as(id, "deta_dt", 11, 3, 0, 1.0, Some("deta_dt"))
            .expect("clean the derivative as in the GUI story");
    }

    let noise = NoiseProfile::Gaussian;
    let filter = FilterKind::SavitzkyGolay;
    let methods = [
        IsoconversionalMethod::OFW,
        IsoconversionalMethod::KAS,
        IsoconversionalMethod::Starink,
        IsoconversionalMethod::FriedmanDifferential,
        IsoconversionalMethod::Vyazovkin,
    ];
    let mut rows = Vec::new();
    for (index, method) in methods.into_iter().enumerate() {
        let output = format!("clean_gate_{index}");
        let (mean_ea, sd_ea, mean_r2, count) =
            run_isoconversional_gui(&model, &ids, method, &output).unwrap_or_else(|error| {
                panic!(
                    "{} failed clean cold-start gate: {error}",
                    method.display_name()
                )
            });
        let relative_error = (mean_ea - E_TRUE).abs() / E_TRUE;
        assert!(
            relative_error <= 0.08,
            "{} failed clean cold-start recovery: mean Ea={mean_ea:.1} J/mol, truth={E_TRUE:.1}, relative error={:.2}%",
            method.display_name(),
            100.0 * relative_error
        );
        let mut row = row_for_result(
            ThermalScenario::NonIsothermal,
            noise,
            seed,
            filter,
            method.display_name(),
            (method != IsoconversionalMethod::Vyazovkin).then_some(mean_r2),
            Some(mean_ea),
            Some(sd_ea),
            None,
            &format!("{count} conversion levels; 8% Ea gate passed"),
        );
        row.noise_seed = format!("noise disabled / {seed}");
        row.filter = "SG 31/3 mass, eta; SG 11/3 rate".into();
        rows.push(row);
    }

    println!(
        "\n=== STORY GATE / CLEAN COLD-START NON-ISOTHERMAL (truth E={E_TRUE:.0} J/mol) ===\n{}",
        Table::new(rows).with(Style::rounded())
    );
}
