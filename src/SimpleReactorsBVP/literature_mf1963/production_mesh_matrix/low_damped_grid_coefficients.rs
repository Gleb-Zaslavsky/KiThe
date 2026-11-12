use super::{
    clustered_input, damped_settings, low_case, mf_discrete_trapezoid_defects,
    passes_mflow_reference_gate, production_profile_errors, try_solve_production,
};
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;

fn compare_coefficients(
    name: &str,
    baseline_method: GridRefinementMethod,
    stricter_method: GridRefinementMethod,
) {
    let case = low_case();
    let (mesh, guess) = clustered_input();
    for (variant, method) in [
        ("baseline", baseline_method),
        ("stricter-gradient-indicator", stricter_method),
    ] {
        let settings = damped_settings(method, 4);
        match try_solve_production(case, mesh.clone(), guess.clone(), settings) {
            Err(error) => {
                println!("MF-LOW grid coefficients {name}/{variant}: SOLVE_FAIL {error}")
            }
            Ok(reactor) => {
                let solution = reactor
                    .solver
                    .get_solution()
                    .expect("successful coefficient run publishes a solution");
                let errors = production_profile_errors(case, &reactor);
                let defects = mf_discrete_trapezoid_defects(case, &reactor);
                let transition_nodes = solution
                    .column(2)
                    .iter()
                    .filter(|value| (0.05..=0.95).contains(*value))
                    .count();
                let gate = passes_mflow_reference_gate(errors, &solution.transpose());
                let score = (errors.0 / 0.5)
                    .max(errors.1 / 2e-3)
                    .max(errors.2 / 2e-3)
                    .max(errors.3 / 2e-8);

                println!(
                    "MF-LOW grid coefficients {name}/{variant}: nodes={} transition_nodes(C0=0.05..0.95)={} dT={:.6e} K dY={:.6e} rel_dq={:.6e} closure={:.3e} defects={defects:?} gate={} score={score:.6e}",
                    solution.nrows(),
                    transition_nodes,
                    errors.0,
                    errors.1,
                    errors.2,
                    errors.3,
                    if gate { "PASS" } else { "FAIL" },
                );
                assert!(solution.iter().all(|value| value.is_finite()));
                assert!(defects.0 < 1e-5 && defects.1 < 1e-5 && defects.2 < 1e-5);
            }
        }
    }
}

#[test]
#[ignore = "MF-LOW TwoPoint refinement-coefficient sensitivity experiment"]
fn mf1963_mflow_damped_twopoint_coefficient_sweep() {
    compare_coefficients(
        "TwoPoint",
        GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0),
        GridRefinementMethod::TwoPoint(0.01, 0.0001, 2.0),
    );
}

#[test]
#[ignore = "MF-LOW Pearson refinement-coefficient sensitivity experiment"]
fn mf1963_mflow_damped_pearson_coefficient_sweep() {
    compare_coefficients(
        "Pearson",
        GridRefinementMethod::Pearson(1e-3, 1.4),
        GridRefinementMethod::Pearson(5e-5, 1.4),
    );
}

#[test]
#[ignore = "MF-LOW Grcar-Smooke refinement-coefficient sensitivity experiment"]
fn mf1963_mflow_damped_grcar_coefficient_sweep() {
    compare_coefficients(
        "Grcar-Smooke",
        GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4),
        GridRefinementMethod::GrcarSmooke(1e-3, 1e-4, 1.4),
    );
}
