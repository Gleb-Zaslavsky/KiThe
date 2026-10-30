use super::SimpleReactorIVP::SimpleReactorTask;
use super::solver_backend::ReactorIvpSolverConfig;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
use std::collections::HashMap;

pub(crate) fn multistep_oxidation_coverage_task() -> SimpleReactorTask {
    let mut task = SimpleReactorTask::new();
    let reactions = vec![
        ReactionData::new_elementary(
            "2CH4+3O2=>2CO+4H2O".to_string(),
            vec![2.5e2, 0.0, 8.0e3],
            None,
        ),
        ReactionData::new_elementary("2CO+O2=>2CO2".to_string(), vec![1.5e2, 0.0, 5.0e3], None),
    ];
    task.kindata
        .set_reaction_data_directly(reactions, None)
        .expect("balanced multistep fixture should parse");
    task.kindata.substances = vec![
        "CH4".into(),
        "O2".into(),
        "CO".into(),
        "H2O".into(),
        "CO2".into(),
    ];
    task.kindata.stecheodata.vec_of_molmasses = Some(vec![16.043, 31.998, 28.010, 18.015, 44.009]);
    // Synthetic rate and heat parameters exercise an intermediate and a
    // consecutive reaction; these are not an experimental mechanism.
    task.thermal_effects = vec![-2.0e5, -2.83e5];
    task.ro = 1350.0;
    task.Cp = 1100.0;
    task.Lambda = 0.22;
    task.m = 0.012;
    task.scaling = ScalingConfig::new(80.0, 0.04, 90.0);
    task.initial_conditions = HashMap::from([
        ("T".into(), 520.0),
        ("q".into(), 0.0),
        ("CH4".into(), 0.25),
        ("O2".into(), 0.75),
        ("CO".into(), 0.0),
        ("H2O".into(), 0.0),
        ("CO2".into(), 0.0),
    ]);
    task.set_solver_backend_config(
        ReactorIvpSolverConfig::default()
            .with_integration_domain(0.0, 0.02)
            .with_max_step(2.0e-4)
            .with_rtol(1.0e-7)
            .with_atol(1.0e-10),
    );
    task
}
