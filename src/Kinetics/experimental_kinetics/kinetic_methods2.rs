use std::fmt::format;

use super::super::solid_state_kinetics_IVP::{IVPSolution, KineticModelIVP, KineticModelNames};
use super::kinetic_methods::integral_isoconversion::IntegralIsoconversionalSolver;
use super::kinetic_methods::{ConversionGridBuilder, ExperimentData, KineticDataView};
use super::one_experiment_dataset::TGADomainError;
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use RustedSciThe::numerical::ODE_api2::{SolverType, UniversalODESolver};
use std::time::Instant;
pub fn solution_to_experiment(sol: IVPSolution, heating_rate: f64, id: &str) -> ExperimentData {
    ExperimentData {
        meta: ExperimentMeta {
            id: id.to_string(),
            heating_rate: Some(heating_rate),
            isothermal_temperature: None,
            comment: None,
        },

        time: sol.time,
        temperature: sol.temperature,
        conversion: sol.conversion,
        conversion_rate: sol.conversion_rate,

        mass: None,
        mass_rate: None,
    }
}

pub fn simulate_tga_dataset(
    model: KineticModelNames,
    e: f64,
    a: f64,
    t0: f64,
    betas: &[f64],
    t_end: f64,
    solvertype: SolverType,
    params: Vec<f64>,
    id: &str,
) -> Result<KineticDataView, TGADomainError> {
    let mut experiments = Vec::with_capacity(betas.len());

    for &beta in betas {
        let mut ivp = KineticModelIVP::new(solvertype.clone());

        let _ = ivp.set_model(model.clone(), params.clone());
        let _ = ivp.set_arrhenius(e, a);
        let _ = ivp.set_heating_program(t0, beta);
        let _ = ivp.set_time(t_end);

        let _ = ivp
            .solve()
            .map_err(|s| TGADomainError::InvalidOperation(format!("error solving ivp {}", s)));

        let sol: IVPSolution = ivp.get_solution().unwrap();

        let exp = solution_to_experiment(sol, beta, id);

        experiments.push(exp);
    }

    Ok(KineticDataView { experiments })
}

#[test]
fn test_conversion_grid_build() {
    use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
    let now = Instant::now();
    let data = simulate_tga_dataset(
        KineticModelNames::F2,
        100_000.0,
        1e6,
        420.0,
        &[0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
        200.0,
        SolverType::Radau(RadauOrder::Order7),
        vec![],
        "",
    )
    .unwrap();
    println!("{:?}", data);
    println!("tga simulated {:?}", now.elapsed());
    let grid = ConversionGridBuilder::new()
        .segments(80)
        .build_nonisothermal(&data)
        .unwrap();

    assert_eq!(grid.temperature.nrows(), 3);
    assert_eq!(grid.eta.len(), 80);

    assert!(grid.temperature[[0, 10]] > 300.0);
}

#[test]
fn test_kas_activation_energy() {
    let true_e = 120_000.0;
    let now = Instant::now();
    let data = simulate_tga_dataset(
        KineticModelNames::F1,
        true_e,
        1e13,
        300.0,
        &[2.0, 5.0, 10.0, 20.0],
        2000.0,
        SolverType::BDF,
        vec![],
        "",
    )
    .unwrap();
    println!("tga simulated {:?}", now.elapsed());
    let grid = ConversionGridBuilder::new()
        .segments(80)
        .build_isothermal(&data)
        .unwrap();

    let solver = IntegralIsoconversionalSolver::kas();
    grid.report();
    let result = solver.solve(&grid).unwrap();
    result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
}

/*
#[test]
fn test_kissinger_method() {

    let true_e = 100_000.0;

    let data = simulate_tga_dataset(
        KineticModel::F1,
        true_e,
        1e13,
        300.0,
        &[2.0, 5.0, 10.0, 20.0],
        2000.0,
    ).unwrap();

    let grid = ConversionGridBuilder::new()
        .segments(200)
        .build(&data)
        .unwrap();

    let res = kissinger_from_grid(&grid).unwrap();

    assert!((res.ea - true_e).abs() < 5000.0);
}

#[test]
fn test_cka_model_detection() {

    let data = simulate_tga_dataset(
        KineticModel::F1,
        120_000.0,
        1e13,
        300.0,
        &[10.0],
        2000.0,
    ).unwrap();

    let solver = CombinedKineticAnalysis::default();

    let res = solver.solve_all(&data).unwrap();

    assert!(res[0].regression.r.abs() > 0.99);
}

#[test]
fn test_conversion_monotonicity() {

    let data = simulate_tga_dataset(
        KineticModel::F1,
        100_000.0,
        1e13,
        300.0,
        &[10.0],
        2000.0,
    ).unwrap();

    let exp = &data.experiments[0];

    assert!(exp.conversion
        .windows(2)
        .all(|w| w[0] <= w[1]));
}

 */
