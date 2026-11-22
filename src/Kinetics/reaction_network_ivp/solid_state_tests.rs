use crate::Kinetics::solid_state_kinetics_IVP::{
    KineticModel, KineticModelIVP, KineticModelNames, SolidStateIvpSolverOptions,
};
use RustedSciThe::numerical::StiffAPI::{StiffFrontend, StiffSolverKind};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::{DMatrix, DVector};

use super::*;

fn stage(
    id: &str,
    reactant: &str,
    product: &str,
    model: KineticModel,
    alpha0: f64,
) -> ReactionSpec {
    ReactionSpec {
        id: id.into(),
        reactants: vec![StoichTerm::new(reactant, 1.0)],
        products: vec![StoichTerm::new(product, 1.0)],
        rate_law: RateLaw::SolidState {
            coefficient: RateCoefficient::Constant(0.5),
            model,
            capacity: 1.0,
            initial_conversion: alpha0,
        },
    }
}

fn problem(reactions: Vec<ReactionSpec>) -> Result<NetworkProblem, NetworkIvpError> {
    NetworkProblem::new(
        vec![
            SpeciesSpec {
                id: "A".into(),
                initial_concentration: 1.0,
            },
            SpeciesSpec {
                id: "B".into(),
                initial_concentration: 0.0,
            },
            SpeciesSpec {
                id: "C".into(),
                initial_concentration: 0.0,
            },
        ],
        reactions,
        TemperatureProgram::Constant(700.0),
        0.0,
        1.0,
    )
}

#[test]
fn f1_stage_matches_closed_form_and_scalar_rate_across_families() {
    let problem = problem(vec![stage("s", "A", "B", KineticModel::F1, 0.1)]).unwrap();
    let assembly = problem.assemble();
    assert_eq!(assembly.variables, ["c0", "c1", "c2", "alpha0"]);
    let rhs = problem
        .evaluate_rhs_with_conversions(0.0, &[1.0, 0.0, 0.0], &[0.1])
        .unwrap();
    assert!((rhs[0] + 0.45).abs() < 1e-12);
    assert!((rhs[1] - 0.45).abs() < 1e-12);
    assert!((rhs[3] - 0.45).abs() < 1e-12);
    for frontend in [StiffFrontend::ExprLegacy, StiffFrontend::AtomViewNative] {
        for family in [
            StiffSolverKind::Bdf,
            StiffSolverKind::Radau,
            StiffSolverKind::BackwardEuler,
            StiffSolverKind::Lsode2,
        ] {
            let options = NetworkIvpOptions {
                family,
                frontend,
                max_step: 0.02,
                ..Default::default()
            };
            let solution = solve_network(&problem, &options).unwrap();
            let alpha = solution.conversions.last().unwrap()[0];
            let expected = 1.0 - 0.9 * (-0.5_f64).exp();
            assert!(
                (alpha - expected).abs() < 2e-3,
                "{family:?}/{frontend:?}: {alpha}"
            );
            assert!((solution.concentrations.last().unwrap()[0] - (1.1 - expected)).abs() < 2e-3);
            assert_eq!(solution.solid_state_reaction_ids, ["s"]);
        }
    }
}

#[test]
fn mixed_noncompeting_network_and_prepared_rebind_preserve_conversion_state() {
    let mut reactions = vec![stage("s", "A", "B", KineticModel::F1, 0.1)];
    reactions.push(ReactionSpec {
        id: "mass".into(),
        reactants: vec![StoichTerm::new("B", 1.0)],
        products: vec![StoichTerm::new("C", 1.0)],
        rate_law: RateLaw::MassAction {
            coefficient: RateCoefficient::Constant(0.2),
        },
    });
    let problem = problem(reactions).unwrap();
    let options = NetworkIvpOptions {
        max_step: 0.02,
        ..Default::default()
    };
    let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
    let first = prepared.solve().unwrap().conversions.last().unwrap()[0];
    assert!(first > 0.1);
    let second = prepared
        .solve_with_initial_concentrations(&[2.0, 0.0, 0.0])
        .unwrap();
    assert!((second.conversions[0][0] - 0.1).abs() < 1e-12);
    let endpoint = second.concentrations.last().unwrap();
    assert!((endpoint.iter().sum::<f64>() - 2.0).abs() < 1e-6);
}

#[test]
fn shared_consumption_and_insufficient_stock_are_typed() {
    let mut reactions = vec![stage("s", "A", "B", KineticModel::F1, 0.1)];
    reactions.push(ReactionSpec {
        id: "mass".into(),
        reactants: vec![StoichTerm::new("A", 1.0)],
        products: vec![StoichTerm::new("C", 1.0)],
        rate_law: RateLaw::MassAction {
            coefficient: RateCoefficient::Constant(0.2),
        },
    });
    assert!(matches!(
        problem(reactions),
        Err(NetworkIvpError::SharedSolidStateResource { .. })
    ));
    let mut too_large = stage("s", "A", "B", KineticModel::F1, 0.1);
    if let RateLaw::SolidState { capacity, .. } = &mut too_large.rate_law {
        *capacity = 2.0;
    }
    assert!(matches!(
        problem(vec![too_large]),
        Err(NetworkIvpError::SolidStateResourceShortfall { .. })
    ));
}

#[test]
fn model_domains_and_negative_rate_parameters_are_rejected_without_clipping() {
    assert!(matches!(
        problem(vec![stage("s", "A", "B", KineticModel::D1, 0.0)]),
        Err(NetworkIvpError::SolidStateDomain { .. })
    ));
    assert!(matches!(
        problem(vec![stage(
            "s",
            "A",
            "B",
            KineticModel::SBtp(1.0, 1.0, -1.0),
            0.5
        )]),
        Err(NetworkIvpError::SolidStateInput { .. })
    ));
    let valid = problem(vec![stage("s", "A", "B", KineticModel::F1, 0.1)]).unwrap();
    assert!(matches!(
        valid.evaluate_rhs_with_conversions(0.0, &[1.0, 0.0, 0.0], &[1.1]),
        Err(NetworkIvpError::SolidStateDomain { .. })
    ));
}

#[test]
fn prepared_solid_state_rate_rebind_and_shortfall_preserve_last_good_result() {
    let problem = problem(vec![stage("s", "A", "B", KineticModel::F1, 0.1)]).unwrap();
    let options = NetworkIvpOptions {
        max_step: 0.02,
        ..Default::default()
    };
    let mut prepared = PreparedNetworkIvp::new(&problem, &options).unwrap();
    let first = prepared.solve().unwrap().conversions.last().unwrap()[0];
    let second = prepared
        .solve_with_rate_coefficients(
            &[1.0, 0.0, 0.0],
            &[RateCoefficientBinding::new(
                "s",
                RateCoefficient::Constant(1.0),
            )],
        )
        .unwrap()
        .conversions
        .last()
        .unwrap()[0];
    assert!(second > first);
    let error = prepared
        .solve_with_initial_concentrations(&[0.1, 0.0, 0.0])
        .unwrap_err();
    assert!(matches!(
        error,
        NetworkIvpError::SolidStateResourceShortfall { .. }
    ));
    assert_eq!(
        prepared
            .last_solution()
            .unwrap()
            .conversions
            .last()
            .unwrap()[0],
        second
    );
}

#[test]
fn invalid_conversion_trajectory_is_not_publishable() {
    let problem = problem(vec![stage("s", "A", "B", KineticModel::F1, 0.1)]).unwrap();
    let axis = DVector::from_vec(vec![0.0, 1.0]);
    let values = DMatrix::from_row_slice(2, 4, &[1.0, 0.0, 0.0, 0.1, 0.9, 0.1, 0.0, 1.1]);
    let error =
        super::solver::validate_candidate(&problem, None, &axis, &values, None).unwrap_err();
    assert!(matches!(error, NetworkIvpError::InvalidConversion { .. }));

    let decreasing = DMatrix::from_row_slice(2, 4, &[0.9, 0.1, 0.0, 0.2, 1.0, 0.0, 0.0, 0.1]);
    let error =
        super::solver::validate_candidate(&problem, None, &axis, &decreasing, None).unwrap_err();
    assert!(matches!(error, NetworkIvpError::InvalidConversion { .. }));
}

#[test]
fn mixed_three_law_network_is_conservative_and_shared_resource_is_rejected_in_either_order() {
    let solid = stage("solid", "A", "B", KineticModel::F1, 0.1);
    let mass = ReactionSpec {
        id: "mass".into(),
        reactants: vec![StoichTerm::new("B", 1.0)],
        products: vec![StoichTerm::new("C", 1.0)],
        rate_law: RateLaw::MassAction {
            coefficient: RateCoefficient::Constant(0.2),
        },
    };
    let custom = ReactionSpec {
        id: "custom".into(),
        reactants: vec![StoichTerm::new("C", 1.0)],
        products: vec![StoichTerm::new("D", 1.0)],
        rate_law: RateLaw::CustomExpr(Expr::Var("c2".into()) * Expr::Const(0.1)),
    };
    let species = vec![
        SpeciesSpec {
            id: "A".into(),
            initial_concentration: 1.0,
        },
        SpeciesSpec {
            id: "B".into(),
            initial_concentration: 0.0,
        },
        SpeciesSpec {
            id: "C".into(),
            initial_concentration: 0.0,
        },
        SpeciesSpec {
            id: "D".into(),
            initial_concentration: 0.0,
        },
    ];
    let make = |reactions| {
        NetworkProblem::new(
            species.clone(),
            reactions,
            TemperatureProgram::Constant(700.0),
            0.0,
            1.0,
        )
    };
    let network = make(vec![solid.clone(), mass, custom.clone()]).unwrap();
    let options = NetworkIvpOptions {
        max_step: 0.02,
        ..Default::default()
    };
    let solution = solve_network(&network, &options).unwrap();
    let last = solution.concentrations.last().unwrap();
    assert!((last.iter().sum::<f64>() - 1.0).abs() < 1e-6);
    assert!(last[1] > 0.0 && last[2] > 0.0 && last[3] > 0.0);

    let competing = ReactionSpec {
        id: "custom_competing".into(),
        reactants: vec![StoichTerm::new("A", 1.0)],
        products: vec![StoichTerm::new("D", 1.0)],
        rate_law: RateLaw::CustomExpr(Expr::Var("c0".into())),
    };
    let second_solid = stage("solid_2", "A", "D", KineticModel::F1, 0.1);
    for reactions in [
        vec![solid.clone(), competing.clone()],
        vec![competing, solid.clone()],
        vec![solid, second_solid],
    ] {
        let error = make(reactions).unwrap_err();
        assert!(matches!(
            error,
            NetworkIvpError::SharedSolidStateResource { .. }
        ));
    }
}

#[test]
fn catalog_model_domains_are_explicit_at_interior_and_endpoints() {
    use super::solid_state::validate_model;

    let models = [
        KineticModel::A2,
        KineticModel::A3,
        KineticModel::A4,
        KineticModel::R2,
        KineticModel::R3,
        KineticModel::D1,
        KineticModel::D2,
        KineticModel::D3,
        KineticModel::D4,
        KineticModel::P2_3,
        KineticModel::P2,
        KineticModel::P3,
        KineticModel::F1,
        KineticModel::F2,
        KineticModel::F3,
        KineticModel::SB(1.0, 1.0, 1.0),
        KineticModel::JMA(2.0),
        KineticModel::Ac(2.0),
        KineticModel::Dec(1.0),
        KineticModel::PTe(1.0, 1.0),
        KineticModel::SBtp(1.0, 1.0, 1.0),
    ];
    for model in &models {
        assert!(validate_model("s", model, 0.5).is_ok(), "{model:?}");
        assert!(validate_model("s", model, -0.1).is_err(), "{model:?}");
        assert!(validate_model("s", model, 1.1).is_err(), "{model:?}");
    }
    for model in [KineticModel::F1, KineticModel::F2, KineticModel::F3] {
        assert!(validate_model("s", &model, 0.0).is_ok());
        assert!(validate_model("s", &model, 1.0).is_ok());
    }
    for model in [KineticModel::R2, KineticModel::R3] {
        assert!(validate_model("s", &model, 0.0).is_ok());
        assert!(validate_model("s", &model, 1.0).is_err());
    }
    for model in [
        KineticModel::P2,
        KineticModel::P3,
        KineticModel::D1,
        KineticModel::P2_3,
    ] {
        assert!(validate_model("s", &model, 0.0).is_err());
        assert!(validate_model("s", &model, 1.0).is_err());
    }
    for model in [
        KineticModel::D1,
        KineticModel::D2,
        KineticModel::D3,
        KineticModel::D4,
        KineticModel::P2_3,
    ] {
        assert!(validate_model("s", &model, 0.0).is_err());
    }
    for model in [
        KineticModel::A2,
        KineticModel::A3,
        KineticModel::A4,
        KineticModel::D2,
        KineticModel::D3,
        KineticModel::D4,
    ] {
        assert!(validate_model("s", &model, 1.0).is_err());
    }
}

#[test]
fn f1_network_matches_existing_scalar_solid_state_route() {
    let alpha0 = 1e-5;
    let mut scalar = KineticModelIVP::new(StiffSolverKind::Bdf);
    scalar.set_problem(1.0, 1.0, 700.0, 5_000.0, 0.5).unwrap();
    scalar.set_model(KineticModelNames::F1, vec![]).unwrap();
    scalar.set_solver_options(SolidStateIvpSolverOptions {
        max_step: 0.01,
        rtol: 1e-8,
        atol: 1e-10,
        ..Default::default()
    });
    scalar.solve().unwrap();
    let scalar_alpha = *scalar.get_solution().unwrap().conversion.last().unwrap();

    let network = NetworkProblem::new(
        vec![
            SpeciesSpec {
                id: "A".into(),
                initial_concentration: 1.0,
            },
            SpeciesSpec {
                id: "B".into(),
                initial_concentration: 0.0,
            },
        ],
        vec![ReactionSpec {
            id: "solid".into(),
            reactants: vec![StoichTerm::new("A", 1.0)],
            products: vec![StoichTerm::new("B", 1.0)],
            rate_law: RateLaw::SolidState {
                coefficient: RateCoefficient::Arrhenius {
                    prefactor: 0.5,
                    activation_energy: 5_000.0,
                },
                model: KineticModel::F1,
                capacity: 1.0,
                initial_conversion: alpha0,
            },
        }],
        TemperatureProgram::Linear {
            initial_kelvin: 700.0,
            kelvin_per_second: 1.0,
        },
        0.0,
        1.0,
    )
    .unwrap();
    let options = NetworkIvpOptions {
        family: StiffSolverKind::Bdf,
        max_step: 0.01,
        ..Default::default()
    };
    let solution = solve_network(&network, &options).unwrap();
    let network_alpha = solution.conversions.last().unwrap()[0];
    assert!(
        (network_alpha - scalar_alpha).abs() < 1e-5,
        "network={network_alpha}, scalar={scalar_alpha}"
    );
    assert!(
        (solution.concentrations.last().unwrap()[0] - (1.0 - (network_alpha - alpha0))).abs()
            < 1e-6
    );
}
