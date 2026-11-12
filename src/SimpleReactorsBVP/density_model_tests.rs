use super::SimpleReactorBVP::DensityModel;
use super::SimpleReactorBVP::{FastElemReact, R_G, SimpleReactorTask};
use super::reactor_BVP_utils::ScalingConfig;
use super::reactor_bvp_test_support::compact_hmx_reactor;
use std::collections::HashMap;

fn prepared(model: DensityModel) -> super::SimpleReactorBVP::SimpleReactorTask {
    let mut reactor = compact_hmx_reactor();
    reactor.set_density_model(model);
    reactor.setup_bvp().expect("density-model BVP setup");
    reactor
}

fn prepared_second_order(model: DensityModel) -> SimpleReactorTask {
    let mut reactor = SimpleReactorTask::new();
    reactor
        .fast_react_set(vec![FastElemReact {
            eq: "A+A=>B".into(),
            A: 2.0,
            n: 0.0,
            E: 0.0,
            Q: 0.0,
        }])
        .expect("second-order reaction definition");
    reactor.kindata.substances = vec!["A".into(), "B".into()];
    reactor.kindata.groups = Some(HashMap::from([
        (
            "A".into(),
            HashMap::from([("C".into(), 2), ("H".into(), 4)]),
        ),
        (
            "B".into(),
            HashMap::from([("C".into(), 4), ("H".into(), 8)]),
        ),
    ]));
    reactor.set_parameters(
        vec![0.0],
        101_325.0,
        500.0,
        1000.0,
        HashMap::from([("T".into(), 500.0), ("A".into(), 0.5), ("B".into(), 0.5)]),
        0.1,
        HashMap::from([("A".into(), 1e-4), ("B".into(), 1e-4)]),
        0.01,
        ScalingConfig::new(350.0, 0.01, 100.0),
    );
    reactor.M = 0.04;
    reactor.set_density_model(model);
    reactor.setup_bvp().expect("second-order BVP setup");
    reactor
}

fn assert_close(actual: f64, expected: f64, tolerance: f64, label: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{label}: actual={actual:.16e}, expected={expected:.16e}, error={:.3e}",
        (actual - expected).abs()
    );
}

#[test]
fn variable_density_is_the_default_and_constant_is_explicit() {
    let mut reactor = compact_hmx_reactor();
    assert_eq!(
        reactor.density_model,
        DensityModel::IdealGasFixedComposition
    );

    reactor.set_density_model(DensityModel::Constant);
    assert_eq!(reactor.density_model, DensityModel::Constant);
}

#[test]
fn variable_density_obeys_reference_law_and_symbolic_derivative() {
    let reactor = prepared(DensityModel::IdealGasFixedComposition);
    let density = reactor
        .density_expression()
        .expect("symbolic density expression");
    let rho = density.lambdify_borrowed_thread_safe(&["Teta"]);
    let derivative = density.diff("Teta");
    let drho = derivative.lambdify_borrowed_thread_safe(&["Teta"]);

    let rho_ref = reactor.ideal_gas_density().expect("reference density");
    let t_scale = reactor.scaling.T_scale;
    let d_t = reactor.scaling.dT;
    let tm = reactor.Tm;
    assert_ne!(d_t, tm, "fixture must distinguish dT from Tm");

    for temperature in [tm * 0.5, tm, tm * 2.0] {
        let teta = (temperature - d_t) / t_scale;
        let expected = rho_ref * tm / temperature;
        assert_close(rho(&[teta]), expected, 2e-13 * expected, "rho(T)");
        assert_close(
            rho(&[teta]) * temperature,
            rho_ref * tm,
            2e-13 * rho_ref * tm,
            "rho(T)*T invariant",
        );
        assert_close(
            drho(&[teta]),
            -rho(&[teta]) * t_scale / temperature,
            2e-12 * expected,
            "d rho / d Teta",
        );

        let h = 1e-5;
        let finite_difference = (rho(&[teta + h]) - rho(&[teta - h])) / (2.0 * h);
        assert_close(
            drho(&[teta]),
            finite_difference,
            2e-8 * drho(&[teta]).abs().max(1.0),
            "central-difference density derivative",
        );
    }
}

#[test]
fn constant_density_matches_the_reference_at_every_temperature() {
    let reactor = prepared(DensityModel::Constant);
    let density = reactor
        .density_expression()
        .expect("constant density expression");
    let rho = density.lambdify_borrowed_thread_safe(&["Teta"]);
    let rho_ref = reactor.ideal_gas_density().expect("reference density");

    for teta in [-0.25, 0.0, 0.75, 2.0, 4.0] {
        assert_close(rho(&[teta]), rho_ref, 0.0, "constant density");
    }
}

#[test]
fn first_order_chemical_rate_scales_with_variable_density() {
    let constant = prepared(DensityModel::Constant);
    let variable = prepared(DensityModel::IdealGasFixedComposition);
    let equation = "HMX=>DecompositionProduct";
    let constant_rate = constant.map_eq_rate.get(equation).expect("constant rate");
    let variable_rate = variable.map_eq_rate.get(equation).expect("variable rate");
    let args = ["Teta", "C0", "C1"];
    let constant_eval = constant_rate.lambdify_borrowed_thread_safe(&args);
    let variable_eval = variable_rate.lambdify_borrowed_thread_safe(&args);

    let temperature = 1800.0;
    let teta = (temperature - variable.scaling.dT) / variable.scaling.T_scale;
    let state = [teta, 0.8, 0.2];
    let ratio = variable_eval(&state) / constant_eval(&state);
    assert_close(
        ratio,
        variable.Tm / temperature,
        2e-12,
        "first-order rate density ratio",
    );
}

#[test]
fn second_order_chemical_rate_has_the_squared_density_ratio() {
    let constant = prepared_second_order(DensityModel::Constant);
    let variable = prepared_second_order(DensityModel::IdealGasFixedComposition);
    let equation = "A+A=>B";
    let constant_rate = constant.map_eq_rate.get(equation).expect("constant rate");
    let variable_rate = variable.map_eq_rate.get(equation).expect("variable rate");
    let args = ["Teta", "C0", "C1"];
    let constant_eval = constant_rate.lambdify_borrowed_thread_safe(&args);
    let variable_eval = variable_rate.lambdify_borrowed_thread_safe(&args);

    let temperature = 750.0;
    let teta = (temperature - variable.scaling.dT) / variable.scaling.T_scale;
    let state = [teta, 0.8, 0.2];
    let ratio = variable_eval(&state) / constant_eval(&state);
    assert_close(
        ratio,
        (variable.Tm / temperature).powi(2),
        2e-12,
        "second-order rate density ratio",
    );
}

#[test]
fn variable_density_arrhenius_rate_has_the_expected_temperature_derivative() {
    let reactor = prepared(DensityModel::IdealGasFixedComposition);
    let equation = "HMX=>DecompositionProduct";
    let rate = reactor.map_eq_rate.get(equation).expect("Arrhenius rate");
    let state_names = ["Teta", "C0", "C1"];
    let rate_eval = rate.lambdify_borrowed_thread_safe(&state_names);
    let derivative = rate.diff("Teta");
    let derivative_eval = derivative.lambdify_borrowed_thread_safe(&state_names);

    let temperature = 1800.0;
    let teta = (temperature - reactor.scaling.dT) / reactor.scaling.T_scale;
    let state = [teta, 0.8, 0.2];
    let rate_value = rate_eval(&state);
    let activation_energy = 5000.0 * 4.184;
    let expected = reactor.scaling.T_scale
        * rate_value
        * (activation_energy / (R_G * temperature.powi(2)) - 1.0 / temperature);
    assert_close(
        derivative_eval(&state),
        expected,
        2e-11 * expected.abs().max(1.0),
        "Arrhenius source derivative with variable density",
    );

    let h = 1e-5;
    let finite_difference = (rate_eval(&[teta + h, state[1], state[2]])
        - rate_eval(&[teta - h, state[1], state[2]]))
        / (2.0 * h);
    assert_close(
        derivative_eval(&state),
        finite_difference,
        2e-7 * derivative_eval(&state).abs().max(1.0),
        "Arrhenius source central difference",
    );
}

#[test]
fn density_selection_does_not_change_transport_snapshots() {
    let constant = prepared(DensityModel::Constant);
    let variable = prepared(DensityModel::IdealGasFixedComposition);

    assert_eq!(constant.D_ro_map, variable.D_ro_map);
    assert_eq!(constant.Pe_D, variable.Pe_D);
    assert_eq!(constant.Pe_q, variable.Pe_q);
    assert_eq!(constant.Cp, variable.Cp);
    assert_eq!(constant.Lambda, variable.Lambda);

    let constant_transport_rhs = &constant.solver.eq_system[2];
    let variable_transport_rhs = &variable.solver.eq_system[2];
    assert_eq!(
        format!("{constant_transport_rhs}"),
        format!("{variable_transport_rhs}"),
        "the species concentration equation must retain the same rho*D closure"
    );
}
