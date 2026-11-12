use super::*;
use RustedSciThe::numerical::BVP_Damp::grid_api::GridRefinementMethod;

#[test]
#[ignore = "paired MF-HIGH direct-RST versus production-facade parity matrix"]
fn mf1963_mfhigh_production_facade_matches_standalone_rst_solver_outputs() {
    let case = super::super::cases()[1];
    let nodes = 301;
    let mesh = (0..nodes)
        .map(|node| node as f64 / (nodes - 1) as f64)
        .collect::<Vec<_>>();
    let guess = super::super::reference_seed(case, nodes);

    let mut equation_probe = super::super::reactor_for(case);
    equation_probe
        .setup_bvp()
        .expect("MF-HIGH production equation setup");
    super::assert_equation_parity(&equation_probe, &mesh, &guess, case);
    let equations = equation_probe.solver.eq_system.clone();

    let routes: Vec<(&str, Box<dyn Fn() -> BvpSolverSettings>)> = vec![
        (
            "BVP_sci/Sparse/tol=1e-6/physical-seed",
            Box::new(|| super::sci_settings_with_tolerance(20_000, 1e-6)),
        ),
        (
            "CorrectionNewton/BandedAuto",
            Box::new(super::correction_newton_settings),
        ),
        (
            "Damped/TwoPoint(0.01,0.01,2.0)/2",
            Box::new(|| super::damped_settings(GridRefinementMethod::TwoPoint(0.01, 0.01, 2.0), 2)),
        ),
        (
            "Damped/Pearson(1e-3,1.4)/2",
            Box::new(|| super::damped_settings(GridRefinementMethod::Pearson(1e-3, 1.4), 2)),
        ),
        (
            "Damped/GrcarSmooke(1e-3,1e-2,1.4)/2",
            Box::new(|| {
                super::damped_settings(GridRefinementMethod::GrcarSmooke(1e-3, 1e-2, 1.4), 2)
            }),
        ),
        (
            "Damped/DoublePoints/4",
            Box::new(|| super::damped_settings(GridRefinementMethod::DoublePoints, 4)),
        ),
        ("Frozen/default", Box::new(|| super::frozen_settings(None))),
        (
            "Frozen/at_high_morm=1e-6",
            Box::new(|| {
                super::frozen_settings(Some(HashMap::from([(
                    "at_high_morm".to_owned(),
                    Some(vec![1e-6]),
                )])))
            }),
        ),
        (
            "Frozen/every_m=1",
            Box::new(|| {
                super::frozen_settings(Some(HashMap::from([(
                    "every_m".to_owned(),
                    Some(vec![1.0]),
                )])))
            }),
        ),
    ];

    for (label, make_settings) in routes {
        let settings = make_settings();
        let (direct_mesh, direct_values) = super::solve_standalone_problem(
            mesh.clone(),
            guess.clone(),
            equations.clone(),
            settings,
        );
        let reactor = super::solve_production(case, mesh.clone(), guess.clone(), make_settings());
        let state_delta =
            super::assert_solution_parity(label, &direct_mesh, &direct_values, &mesh, &reactor);
        let direct_errors =
            super::standalone_reference_errors(case, &direct_mesh, &direct_values, reactor.L);
        let facade_errors = super::production_profile_errors(case, &reactor);
        let direct_pass = super::passes_mflow_reference_gate(direct_errors, &direct_values);
        let facade_values = reactor
            .solver
            .get_solution()
            .expect("published production solution");
        let facade_pass =
            super::passes_mflow_reference_gate(facade_errors, &facade_values.transpose());

        println!(
            "MF-HIGH parity {label}: direct_nodes={} facade_nodes={} max_state_scaled={state_delta:?} direct_gate={direct_pass} facade_gate={facade_pass} direct(dT={:.6e},dY={:.6e},rel_dq={:.6e}) facade(dT={:.6e},dY={:.6e},rel_dq={:.6e})",
            direct_mesh.len(),
            reactor.solver.x_mesh.as_ref().unwrap().len(),
            direct_errors.0,
            direct_errors.1,
            direct_errors.2,
            facade_errors.0,
            facade_errors.1,
            facade_errors.2,
        );
        assert!(
            direct_pass,
            "standalone RST missed MF-HIGH gate for {label}"
        );
        assert!(
            facade_pass,
            "production facade missed MF-HIGH gate for {label}"
        );
        assert_eq!(
            facade_pass, direct_pass,
            "gate classification differs for {label}"
        );
        assert!(
            state_delta.iter().copied().fold(0.0_f64, f64::max) <= 2e-10,
            "production facade changed the direct-RST profile for {label}: {state_delta:?}"
        );
    }

    // The second production BVP_sci seed is also passed through the same direct
    // RST/facade comparison, independently of the solver matrix above.
    let smooth_guess = super::super::smooth_seed(case, nodes);
    let settings = super::sci_settings_with_tolerance(20_000, 1e-6);
    let (direct_mesh, direct_values) =
        super::solve_standalone_problem(mesh.clone(), smooth_guess.clone(), equations, settings);
    let reactor = super::solve_production(
        case,
        mesh.clone(),
        smooth_guess,
        super::sci_settings_with_tolerance(20_000, 1e-6),
    );
    let state_delta = super::assert_solution_parity(
        "BVP_sci/Sparse/tol=1e-6/smooth-seed",
        &direct_mesh,
        &direct_values,
        &mesh,
        &reactor,
    );
    assert!(
        state_delta.iter().copied().fold(0.0_f64, f64::max) <= 2e-10,
        "production facade changed the smooth-seed BVP_sci profile: {state_delta:?}"
    );
}
