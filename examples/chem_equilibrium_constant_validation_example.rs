//! Guide: independent equilibrium-constant validation.
//!
//! This example shows the small-system validator path used as a second opinion
//! for `ln(Q) - ln(K)` checks.

use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_problem::{
    EquilibriumConstantActivityModel, EquilibriumConstantProblem,
};
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_solver::EquilibriumConstantSolver;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_constant_validation::{
    EquilibriumConstantValidationTolerances, validate_equilibrium_constants,
};
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_reaction_basis::{
    ReactionBasisTolerances, ValidatedReactionBasis,
};
use nalgebra::DMatrix;
use std::rc::Rc;

fn main() {
    let basis = ValidatedReactionBasis::new(
        vec!["A2".to_string(), "A".to_string()],
        &DMatrix::from_row_slice(2, 1, &[2.0, 1.0]),
        DMatrix::from_column_slice(2, 1, &[-1.0, 2.0]),
        1,
        ReactionBasisTolerances::default(),
    )
    .expect("failed to build reaction basis");

    let problem = EquilibriumConstantProblem::new(
        basis,
        vec![1.0, 0.0],
        vec![Rc::new(|_| -50_000.0), Rc::new(|_| 0.0)],
        EquilibriumConditions::new(1500.0, 2.0 * 101_325.0, 101_325.0)
            .expect("invalid thermodynamic conditions"),
        EquilibriumConstantActivityModel::IdealGas,
    )
    .expect("failed to build validation problem");

    let solution = EquilibriumConstantSolver::default()
        .solve(&problem)
        .expect("independent validation solve failed");
    let validation = validate_equilibrium_constants(
        &problem,
        &solution.moles,
        EquilibriumConstantValidationTolerances::default(),
    )
    .expect("independent validation failed");

    println!("extent = {:.6}", solution.extent);
    println!("accepted = {}", validation.accepted);
    println!(
        "max |ln(Q)-ln(K)| = {:.3e}",
        validation.max_abs_log_residual
    );
}
