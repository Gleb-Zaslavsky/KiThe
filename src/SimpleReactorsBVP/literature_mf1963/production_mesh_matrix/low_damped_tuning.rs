use super::{
    clustered_input, damped_settings_with_tolerance, low_case, mf_discrete_trapezoid_defects,
    passes_mflow_reference_gate, production_profile_errors, try_solve_production,
};
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;

#[test]
#[ignore = "MF-LOW Damped tolerance and refinement-threshold sensitivity sweep"]
fn mf1963_mflow_damped_tolerance_and_grid_threshold_sweep() {
    let case = low_case();
    let (mesh, guess) = clustered_input();
    let methods = [
        (
            "TwoPoint",
            GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0),
            GridRefinementMethod::TwoPoint(0.005, 0.005, 1.4),
        ),
        (
            "Pearson",
            GridRefinementMethod::Pearson(1e-3, 1.4),
            GridRefinementMethod::Pearson(1e-4, 1.4),
        ),
        (
            "Grcar-Smooke",
            GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4),
            GridRefinementMethod::GrcarSmooke(1e-4, 1e-3, 1.4),
        ),
    ];
    let mut best_gate_score = f64::INFINITY;

    for (name, baseline_method, tuned_method) in methods {
        let mut baseline_score = None;
        let mut tuned_score = None;
        for (variant, method, tolerance) in [
            ("baseline-atol-rtol", baseline_method.clone(), 2e-7),
            ("tight-atol-rtol", baseline_method, 1e-9),
            ("tuned-indicator", tuned_method, 2e-7),
        ] {
            let settings = damped_settings_with_tolerance(method, 4, tolerance);
            match try_solve_production(case, mesh.clone(), guess.clone(), settings) {
                Err(error) => {
                    println!("MF-LOW Damped tuning {name}/{variant}: SOLVE_FAIL {error}");
                }
                Ok(reactor) => {
                    let solution = reactor
                        .solver
                        .get_solution()
                        .expect("successful tuning run publishes a solution");
                    let errors = production_profile_errors(case, &reactor);
                    let defects = mf_discrete_trapezoid_defects(case, &reactor);
                    let gate = passes_mflow_reference_gate(errors, &solution.transpose());
                    let score = (errors.0 / 0.5)
                        .max(errors.1 / 2e-3)
                        .max(errors.2 / 2e-3)
                        .max(errors.3 / 2e-8);
                    best_gate_score = best_gate_score.min(score);
                    match variant {
                        "baseline-atol-rtol" => baseline_score = Some(score),
                        "tuned-indicator" => tuned_score = Some(score),
                        _ => {}
                    }
                    println!(
                        "MF-LOW Damped tuning {name}/{variant}: tol={tolerance:.1e} nodes={} dT={:.6e} K dY={:.6e} rel_dq={:.6e} closure={:.3e} defects={defects:?} gate={} score={score:.6e}",
                        solution.nrows(),
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
        assert!(
            tuned_score.expect("tuned indicator result") < baseline_score.expect("baseline result"),
            "the tested stricter {name} indicator should improve the normalized reference error"
        );
    }

    println!("MF-LOW Damped tuning best normalized reference-gate score={best_gate_score:.6e}");
}
