//! Solve a user-defined homogeneous reaction network with a prescribed
//! temperature program. Run with `cargo run --example reaction_network_ivp`.

use KiThe::Kinetics::reaction_network_ivp::{
    NetworkIvpOptions, NetworkProblem, RateCoefficient, RateLaw, ReactionSpec, SpeciesSpec,
    StoichTerm, TemperatureProgram, solve_network,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let problem = NetworkProblem::new(
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
        vec![
            ReactionSpec {
                id: "A_to_2B".into(),
                reactants: vec![StoichTerm::new("A", 1.0)],
                products: vec![StoichTerm::new("B", 2.0)],
                rate_law: RateLaw::MassAction {
                    coefficient: RateCoefficient::Arrhenius {
                        prefactor: 0.8,
                        activation_energy: 18_000.0,
                    },
                },
            },
            ReactionSpec {
                id: "B_to_C".into(),
                reactants: vec![StoichTerm::new("B", 1.0)],
                products: vec![StoichTerm::new("C", 1.0)],
                rate_law: RateLaw::CustomExpr(
                    RustedSciThe::symbolic::symbolic_engine::Expr::Const(0.15)
                        * RustedSciThe::symbolic::symbolic_engine::Expr::Var("c1".into()),
                ),
            },
        ],
        TemperatureProgram::Linear {
            initial_kelvin: 500.0,
            kelvin_per_second: 20.0,
        },
        0.0,
        10.0,
    )?;

    let solution = solve_network(&problem, &NetworkIvpOptions::default())?;
    let final_state = solution
        .concentrations
        .last()
        .expect("solver returned samples");
    println!(
        "t = {:.3} s, T = {:.2} K",
        solution.time_seconds.last().unwrap(),
        solution.temperature_kelvin.last().unwrap()
    );
    for (species, concentration) in solution.species_ids.iter().zip(final_state) {
        println!("{species}: {concentration:.8e} mol/m^3");
    }
    Ok(())
}
