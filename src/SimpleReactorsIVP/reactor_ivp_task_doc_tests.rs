//! Physical file-to-solution stories for checked-in reactor IVP task documents.

use crate::SimpleReactorsIVP::SimpleReactorIVP::SimpleReactorTask;
use crate::gui::document_lifecycle::GuiFileOperationResult;
use crate::gui::reactor_ivp_gui::ReactorIvpApp;
use std::path::PathBuf;

fn task_doc(relative_path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("task_docs/reactors/ivp")
        .join(relative_path)
}

#[test]
fn checked_in_physical_ivp_documents_load_and_solve() {
    for (file, expected_species, expected_reactions, requires_evolution) in [
        ("h2_dissociation.task", 2, 1, false),
        ("two_step_oxidation.task", 5, 2, false),
        ("h2_dissociation_with_heat_release.task", 2, 1, true),
    ] {
        let mut app = ReactorIvpApp::new();
        let path = task_doc(file);
        assert!(matches!(
            app.load_document_from_path(path.clone()),
            GuiFileOperationResult::Loaded { .. }
        ));
        assert_eq!(app.config.substances.len(), expected_species, "{file}");
        assert_eq!(app.config.reactions.len(), expected_reactions, "{file}");

        let mut task: SimpleReactorTask = app.task.clone();
        app.config
            .apply_to_task(&mut task)
            .unwrap_or_else(|error| panic!("{file}: applying the loaded document failed: {error}"));
        let snapshot = task
            .solve()
            .unwrap_or_else(|error| panic!("{file}: physical IVP solve failed: {error}"));
        assert!(
            snapshot.axis.iter().all(|value| value.is_finite()),
            "{file}"
        );
        assert!(
            snapshot.values.iter().all(|value| value.is_finite()),
            "{file}"
        );

        let temperature_index = snapshot
            .variable_names
            .iter()
            .position(|variable| variable == "Teta")
            .expect("physical IVP state should include scaled temperature");
        let temperatures = snapshot
            .values
            .column(temperature_index)
            .map(|teta| teta * task.scaling.T_scale + task.scaling.dT);
        assert!(
            temperatures
                .iter()
                .all(|temperature| (0.0..=3000.0).contains(temperature)),
            "{file}: physical temperature outside [0, 3000] K: {:.3}..{:.3} K",
            temperatures.iter().copied().fold(f64::INFINITY, f64::min),
            temperatures
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
        );
        if requires_evolution {
            let temperature_span = temperatures
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                - temperatures.iter().copied().fold(f64::INFINITY, f64::min);
            assert!(
                temperature_span > 1e-4,
                "{file}: expected measurable thermal evolution, got {temperature_span} K"
            );
            let h2_species = task
                .kindata
                .substances
                .iter()
                .position(|name| name == "H2")
                .expect("heat-release fixture must contain H2");
            let h2_column = snapshot
                .variable_names
                .iter()
                .position(|name| name == &format!("C{h2_species}"))
                .expect("H2 mass-fraction state should be published");
            let h2_initial = task.initial_conditions["H2"];
            let h2_final = snapshot.values[(snapshot.values.nrows() - 1, h2_column)];
            assert!((h2_final - h2_initial).abs() > 1e-6);
            println!(
                "IVP task-doc {file}: dT={temperature_span:.6} K, H2={h2_initial:.8}->{h2_final:.8}"
            );
            let balances = task
                .latest_conservation_report()
                .expect("heat-release story should publish conservation diagnostics");
            assert!(
                balances.energy_balance_error_rel.is_finite()
                    && balances.energy_balance_error_abs.is_finite(),
                "{file}: energy diagnostics must remain finite"
            );
            assert!(
                balances.energy_balance_error_abs.abs() <= 1.0
                    && balances.energy_balance_error_rel.abs() <= 0.01,
                "{file}: energy balance should close within 1 J and 0.01%, got abs={} rel={} %",
                balances.energy_balance_error_abs,
                balances.energy_balance_error_rel
            );
            assert!(
                balances
                    .atomic_mass_balance_error
                    .iter()
                    .all(|(_, error)| error.abs() <= 5e-2),
                "{file}: elemental balance error exceeded the established envelope"
            );
            assert!(
                balances
                    .sum_of_mass_fractions
                    .iter()
                    .all(|(_, sum)| (sum - 1.0).abs() <= 5e-2),
                "{file}: mass-fraction sum left the established envelope"
            );
            println!(
                "IVP task-doc {file}: energy residual abs={:.6e}, rel={:.6e}% (acceptance)",
                balances.energy_balance_error_abs, balances.energy_balance_error_rel
            );
            let q_index = snapshot
                .variable_names
                .iter()
                .position(|name| name == "q")
                .expect("snapshot should include heat flux");
            let temperature_index = snapshot
                .variable_names
                .iter()
                .position(|name| name == "Teta")
                .expect("snapshot should include scaled temperature");
            let x = snapshot.axis.iter().map(|z| z * task.L).collect::<Vec<_>>();
            let q_physical = snapshot
                .values
                .column(q_index)
                .iter()
                .map(|q| q * task.scaling.T_scale / task.L)
                .collect::<Vec<_>>();
            let variable_names = snapshot
                .variable_names
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>();
            let heat_release = snapshot
                .values
                .row_iter()
                .map(|state| {
                    let state = state.iter().copied().collect::<Vec<_>>();
                    task.heat_release.eval_expression(&variable_names, &state)
                })
                .collect::<Vec<_>>();
            let trapezoid = |values: &[f64]| {
                x.windows(2)
                    .zip(values.windows(2))
                    .map(|(x_pair, y_pair)| (x_pair[1] - x_pair[0]) * (y_pair[0] + y_pair[1]) / 2.0)
                    .sum::<f64>()
            };
            let integrated_heat = trapezoid(&heat_release);
            let integrated_q = trapezoid(&q_physical);
            let delta_q = q_physical.last().unwrap() - q_physical[0];
            let delta_temperature = snapshot.values
                [(snapshot.values.nrows() - 1, temperature_index)]
                * task.scaling.T_scale
                - snapshot.values[(0, temperature_index)] * task.scaling.T_scale;
            let sensible_heat = task.m * task.Cp * delta_temperature;
            let independent_energy_residual = -delta_q - integrated_heat + sensible_heat;
            let q_temperature_residual =
                sensible_heat - task.m * task.Cp / task.Lambda * integrated_q;
            assert!(
                independent_energy_residual.abs() <= 1.0,
                "{file}: independent trapezoidal energy residual exceeded 1 J: {independent_energy_residual}"
            );
            println!(
                "IVP task-doc {file}: nodes={} nonuniform_intervals={} dQ={:.6e} integral_H={:.6e} sensible={:.6e} trapezoid_residual={:.6e} q_temperature_residual={:.6e}",
                x.len(),
                x.windows(2)
                    .filter(|pair| ((pair[1] - pair[0]) - (x[1] - x[0])).abs() > 1e-12)
                    .count(),
                delta_q,
                integrated_heat,
                sensible_heat,
                independent_energy_residual,
                q_temperature_residual
            );
        }

        for (index, variable) in snapshot.variable_names.iter().enumerate() {
            if variable == "Teta" || variable == "q" {
                continue;
            }
            assert!(
                snapshot
                    .values
                    .column(index)
                    .iter()
                    .all(|fraction| (0.0..=1.0).contains(fraction)),
                "{file}: invalid mass fraction for {variable}"
            );
            if requires_evolution {
                let species_name = task
                    .kindata
                    .substances
                    .get(index - 2)
                    .expect("species state columns follow Teta and q");
                let initial = task.initial_conditions[species_name];
                let final_value = snapshot.values[(snapshot.values.nrows() - 1, index)];
                assert!(
                    (final_value - initial).abs() > 1e-6,
                    "{file}: {species_name} did not measurably evolve"
                );
            }
        }
        println!(
            "IVP task-doc {file}: solved {} species / {} reactions, T={:.2}..{:.2} K, points={}",
            expected_species,
            expected_reactions,
            temperatures.iter().copied().fold(f64::INFINITY, f64::min),
            temperatures
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max),
            snapshot.time_points()
        );
    }
}
