//! I1-I2 regression tests for the independent pure-phase P,H validator.
//!
//! These fixtures are intentionally synthetic.  A passing test proves the
//! second mathematical route and its scalar contracts, not production phase
//! lifecycle or independent database provenance.

use std::rc::Rc;
use std::sync::Arc;

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EnthalpyScale, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{
    EquilibriumSolverSettings, GibbsFn, Phase, Solvers,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, SolveError,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_formulation::PreparedPhFormulation;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_monolithic::PreparedMonolithicPhRunner;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_options::PhMonolithicOptions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::MolarThermoFunction;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::{
    ResolvedThermochemistry, ThermochemistryProvenance,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::{
    EquilibriumConditions, EquilibriumProblem, LogMolesInitialGuess, PreparedEquilibriumProblem,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::SolverPolicy;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    MOLAR_GAS_CONSTANT, PurePhaseBoundaryElementComposition, PurePhaseBoundaryStructuralTolerances,
};
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhCanonicalEvidence, PurePhasePhConditions, PurePhasePhCrossValidationTolerances,
    PurePhasePhProblem, PurePhasePhSolverSettings, PurePhasePhValidator,
    compare_pure_phase_ph_validation,
};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};

const T_STAR: f64 = 800.0;
const XI_STAR: f64 = 0.25;
// Keep the van't Hoff shift small enough that the inner chemical root stays
// inside the physical extent interval across the complete 650..950 K bracket.
const SYNTHETIC_REACTION_ENTHALPY: f64 = 2_000.0;

fn function<F>(callback: F) -> MolarThermoFunction
where
    F: Fn(f64) -> Result<f64, ReactionExtentError> + Send + Sync + 'static,
{
    Arc::new(callback)
}

fn synthetic_ln_q(extent: f64) -> f64 {
    let a = 2.0 - 2.0 * extent;
    let b = 2.0 - extent;
    let total = a + b;
    -2.0 * (a / total).ln() - (b / total).ln()
}

fn synthetic_ln_k(temperature: f64) -> f64 {
    synthetic_ln_q(XI_STAR)
        + SYNTHETIC_REACTION_ENTHALPY / MOLAR_GAS_CONSTANT * (1.0 / T_STAR - 1.0 / temperature)
}

/// I1/I2 fixture: `2 A(g) + B(g) <=> S(condensed)`.
///
/// `T_star` and the physical candidate amount are declared first. The target
/// enthalpy is then computed directly from this fixture, so the expected root
/// does not come from the validator under test.
fn synthetic_problem(reaction_scale: f64, target_enthalpy_offset: f64) -> PurePhasePhProblem {
    synthetic_problem_with_initial_candidate(reaction_scale, target_enthalpy_offset, 0.0)
}

fn synthetic_problem_with_initial_candidate(
    reaction_scale: f64,
    target_enthalpy_offset: f64,
    initial_candidate_moles: f64,
) -> PurePhasePhProblem {
    synthetic_problem_with_conditions(
        reaction_scale,
        target_enthalpy_offset,
        initial_candidate_moles,
        1.0,
        101_325.0,
        101_325.0,
    )
}

/// Rebuilds the synthetic physical family at a different inventory or pressure
/// reference without mutating a previous accepted P,H case.
///
/// `inventory_scale` scales every extensive quantity, including the target
/// enthalpy. Scaling pressure and reference pressure together preserves the
/// ideal-gas activity ratio and must therefore preserve the intensive state.
fn synthetic_problem_with_conditions(
    reaction_scale: f64,
    target_enthalpy_offset: f64,
    initial_candidate_moles: f64,
    inventory_scale: f64,
    pressure: f64,
    reference_pressure: f64,
) -> PurePhasePhProblem {
    assert!(
        inventory_scale.is_finite() && inventory_scale > 0.0,
        "synthetic inventory scale must be finite and positive"
    );
    let physical_candidate = XI_STAR * inventory_scale;
    let physical_progress = physical_candidate - initial_candidate_moles;
    assert!(
        physical_progress > 0.0,
        "fixture requires a positive phase-forming extent"
    );
    // The synthetic `ln(K)` obeys van't Hoff for a constant reaction
    // enthalpy. Together with `G_s^0 = -R*T*ln(K)` and `H_s^0 = Delta H`,
    // this makes the independent route thermodynamically consistent with the
    // canonical P,H formulation and its analytical temperature Jacobian.
    let candidate_gibbs = function(move |temperature| {
        Ok(-MOLAR_GAS_CONSTANT * temperature * synthetic_ln_k(temperature))
    });

    // Additive synthetic enthalpy law. It is deliberately supplied as its own
    // capability: this validator is testing the P,H equations, not assuming
    // a particular polynomial representation or a hidden production model.
    // The gas reference enthalpies are zero; the condensed candidate carries
    // the constant reaction enthalpy used by `synthetic_ln_k` above.
    let gas_enthalpies = vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))];
    let candidate_enthalpy = function(|_| Ok(SYNTHETIC_REACTION_ENTHALPY));
    let expected_h = physical_candidate * SYNTHETIC_REACTION_ENTHALPY;
    let conditions = PurePhasePhConditions::new(
        pressure,
        reference_pressure,
        expected_h + target_enthalpy_offset,
        650.0,
        950.0,
    )
    .unwrap();

    PurePhasePhProblem::new(
        vec!["A(g)".into(), "B(g)".into()],
        vec![
            1.5 * inventory_scale + 2.0 * physical_progress,
            1.75 * inventory_scale + physical_progress,
        ],
        vec![-2.0 * reaction_scale, -reaction_scale],
        initial_candidate_moles,
        reaction_scale,
        "S(cond)",
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        candidate_gibbs,
        gas_enthalpies,
        candidate_enthalpy,
        conditions,
    )
    .unwrap()
    .with_element_composition(
        PurePhaseBoundaryElementComposition::new(
            vec!["A".into(), "B".into()],
            DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]),
            vec![2.0, 1.0],
        )
        .unwrap(),
        PurePhaseBoundaryStructuralTolerances::default(),
    )
    .unwrap()
}

/// The same physical A/B/S family with the gas components deliberately
/// permuted. The identity-bearing rows of the independent element matrix are
/// permuted with the component data, so the test cannot pass by comparing two
/// same-length anonymous vectors.
fn permuted_synthetic_problem() -> PurePhasePhProblem {
    let candidate_gibbs =
        function(|temperature| Ok(-MOLAR_GAS_CONSTANT * temperature * synthetic_ln_k(temperature)));
    let expected_h = XI_STAR * SYNTHETIC_REACTION_ENTHALPY;
    PurePhasePhProblem::new(
        vec!["B(g)".into(), "A(g)".into()],
        vec![1.75 + XI_STAR, 1.5 + 2.0 * XI_STAR],
        vec![-1.0, -2.0],
        0.0,
        1.0,
        "S(cond)",
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        candidate_gibbs,
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        function(|_| Ok(SYNTHETIC_REACTION_ENTHALPY)),
        PurePhasePhConditions::new(101_325.0, 101_325.0, expected_h, 650.0, 950.0).unwrap(),
    )
    .unwrap()
    .with_element_composition(
        PurePhaseBoundaryElementComposition::new(
            vec!["A".into(), "B".into()],
            // Rows are B, A; columns remain A, B.
            DMatrix::from_row_slice(2, 2, &[0.0, 1.0, 1.0, 0.0]),
            vec![2.0, 1.0],
        )
        .unwrap(),
        PurePhaseBoundaryStructuralTolerances::default(),
    )
    .unwrap()
}

/// I3 fixed-topology bridge only. This constructs the same synthetic physical
/// family through the canonical coupled `[log-moles, temperature]` runner;
/// there is no active-set lifecycle or `solve_resolved_ph` facade involved.
fn canonical_fixed_topology_evidence(
    independent_problem: &PurePhasePhProblem,
) -> PurePhasePhCanonicalEvidence {
    let initial_moles = vec![1.8, 1.9, 0.1];
    let conditions = EquilibriumConditions::new(T_STAR, 101_325.0, 101_325.0).unwrap();
    let candidate_gibbs: GibbsFn =
        Rc::new(|temperature| -MOLAR_GAS_CONSTANT * temperature * synthetic_ln_k(temperature));
    let canonical_problem = EquilibriumProblem::new(
        vec![
            "A(g)".to_string(),
            "B(g)".to_string(),
            "S(cond)".to_string(),
        ],
        initial_moles.clone(),
        LogMolesInitialGuess::from_moles(&initial_moles, 1e-12).unwrap(),
        DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 0.0, 1.0, 2.0, 1.0]),
        vec![Rc::new(|_| 0.0), Rc::new(|_| 0.0), candidate_gibbs],
        vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![2],
            },
        ],
        conditions,
    )
    .unwrap();
    let prepared = PreparedEquilibriumProblem::new(canonical_problem).unwrap();
    let bounds = TemperatureBounds::new(650.0, 950.0).unwrap();
    let phase_gas = PhaseId::new(Some("gas".into()));
    let phase_condensed = PhaseId::new(Some("condensed".into()));
    let thermochemistry = ResolvedThermochemistry::from_functions(
        vec![
            ThermochemistryProvenance::new(
                PhaseComponentId::new(phase_gas.clone(), "A(g)"),
                "synthetic",
                "A",
                "gas",
            ),
            ThermochemistryProvenance::new(
                PhaseComponentId::new(phase_gas, "B(g)"),
                "synthetic",
                "B",
                "gas",
            ),
            ThermochemistryProvenance::new(
                PhaseComponentId::new(phase_condensed, "S(cond)"),
                "synthetic",
                "S",
                "condensed",
            ),
        ],
        bounds,
        vec![
            function(|_| Ok(0.0)),
            function(|_| Ok(0.0)),
            function(|temperature| {
                Ok(-MOLAR_GAS_CONSTANT * temperature * synthetic_ln_k(temperature))
            }),
        ],
        vec![
            function(|_| Ok(0.0)),
            function(|_| Ok(0.0)),
            function(|_| Ok(SYNTHETIC_REACTION_ENTHALPY)),
        ],
        vec![
            Some(function(|_| Ok(0.0))),
            Some(function(|_| Ok(0.0))),
            Some(function(|_| Ok(0.0))),
        ],
    )
    .unwrap();
    let initial_enthalpies = thermochemistry.evaluate_enthalpy(T_STAR).unwrap();
    let target = independent_problem.conditions().target_enthalpy();
    let scale =
        EnthalpyScale::from_magnitudes(target, &initial_moles, &initial_enthalpies).unwrap();
    let formulation =
        PreparedPhFormulation::new(prepared, thermochemistry, bounds, target, scale).unwrap();
    let mut settings = EquilibriumSolverSettings::default();
    settings.solver = Solvers::LM;
    settings.solver_policy = Some(SolverPolicy::legacy_default(Solvers::LM));
    settings.solver_params.tol = 1e-11;
    settings.solver_params.max_iter = 300;
    let runner =
        PreparedMonolithicPhRunner::new(formulation, settings, PhMonolithicOptions::default())
            .unwrap();
    let outcome = runner.solve_from_temperature_seed(T_STAR).unwrap();
    let snapshot = outcome.snapshot;
    let chemical_log_residual = independent_problem
        .chemical_log_residual_for_gas_moles(&snapshot.moles[..2], snapshot.temperature)
        .unwrap();
    let max_abs_element_balance = [
        (snapshot.moles[0] + 2.0 * snapshot.moles[2] - 2.0).abs(),
        (snapshot.moles[1] + snapshot.moles[2] - 2.0).abs(),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    PurePhasePhCanonicalEvidence::new(
        independent_problem.case_identity(),
        snapshot.temperature,
        snapshot.moles[..2].to_vec(),
        snapshot.moles[2],
        snapshot.total_enthalpy,
        Some(chemical_log_residual),
        Some(max_abs_element_balance),
    )
    .unwrap()
}

#[test]
fn i1_i2_recovers_constructed_temperature_extent_and_conservation() {
    // I1/I2: independent structural validation and constructed scalar root.
    let problem = synthetic_problem(1.0, 0.0);
    let structure = problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .expect("I1 fixture must be structurally one strict phase-forming reaction");
    assert_eq!(structure.full_reaction_dimension, 1);
    assert_eq!(structure.gas_only_reaction_dimension, 0);
    assert_eq!(problem.gas_species(), ["A(g)", "B(g)"]);
    assert_eq!(problem.candidate_name(), "S(cond)");
    assert_eq!(problem.gas_stoichiometry(), [-2.0, -1.0]);
    assert_eq!(problem.candidate_stoichiometry(), 1.0);
    assert_eq!(problem.conditions().pressure(), 101_325.0);
    assert_eq!(problem.conditions().reference_pressure(), 101_325.0);
    assert_eq!(problem.conditions().lower_temperature(), 650.0);
    assert_eq!(problem.conditions().upper_temperature(), 950.0);
    assert!(problem.conditions().target_enthalpy() > 0.0);

    let result = PurePhasePhValidator::default()
        .solve(&problem)
        .expect("I2 nested scalar validator must recover the constructed state");
    assert!((result.temperature - T_STAR).abs() <= 1e-6);
    assert!((result.extent - XI_STAR).abs() <= 1e-9);
    assert!((result.gas_moles[0] - 1.5).abs() <= 1e-9);
    assert!((result.gas_moles[1] - 1.75).abs() <= 1e-9);
    assert!((result.candidate_moles - XI_STAR).abs() <= 1e-9);
    assert!(result.chemical_log_residual.abs() <= 1e-10);
    assert!(result.enthalpy_residual.abs() <= 1e-7);
    assert!(result.total_inner_solves >= 3);
}

#[test]
fn i1_i2_coordinate_scaling_preserves_physical_state() {
    // I2: reaction-coordinate rescaling must not alter the physical state.
    let validator = PurePhasePhValidator::default();
    let reference = validator.solve(&synthetic_problem(1.0, 0.0)).unwrap();
    for scale in [0.5, 1.0, 2.0] {
        let result = validator.solve(&synthetic_problem(scale, 0.0)).unwrap();
        assert!((result.temperature - reference.temperature).abs() <= 1e-6);
        assert!((result.extent * scale - reference.extent).abs() <= 1e-9);
        assert!((result.candidate_moles - reference.candidate_moles).abs() <= 1e-9);
        assert!((result.total_enthalpy - reference.total_enthalpy).abs() <= 1e-7);
        for (actual, expected) in result.gas_moles.iter().zip(&reference.gas_moles) {
            assert!((actual - expected).abs() <= 1e-9);
        }
    }
}

#[test]
fn p10_1_ph_inventory_target_enthalpy_scaling_preserves_intensive_state() {
    let validator = PurePhasePhValidator::default();
    let reference = validator.solve(&synthetic_problem(1.0, 0.0)).unwrap();

    for inventory_scale in [0.25, 1.0, 8.0] {
        let result = validator
            .solve(&synthetic_problem_with_conditions(
                1.0,
                0.0,
                0.0,
                inventory_scale,
                101_325.0,
                101_325.0,
            ))
            .unwrap();
        assert!((result.temperature - reference.temperature).abs() <= 1e-6);
        assert!(
            (result.candidate_moles - inventory_scale * reference.candidate_moles).abs() <= 1e-8
        );
        assert!((result.total_enthalpy - inventory_scale * reference.total_enthalpy).abs() <= 1e-6);
        assert!(result.chemical_log_residual.abs() <= 1e-10);
        assert!(result.enthalpy_residual.abs() <= 1e-6);
        for (actual, expected) in result.gas_moles.iter().zip(&reference.gas_moles) {
            assert!(
                (actual - inventory_scale * expected).abs() <= 1e-8,
                "inventory_scale={inventory_scale:e}, actual={actual:e}, expected={expected:e}"
            );
        }
    }
}

#[test]
fn p10_1_ph_joint_pressure_reference_scaling_preserves_solution() {
    let validator = PurePhasePhValidator::default();
    let reference = validator.solve(&synthetic_problem(1.0, 0.0)).unwrap();

    for pressure_scale in [0.1, 1.0, 10.0] {
        let result = validator
            .solve(&synthetic_problem_with_conditions(
                1.0,
                0.0,
                0.0,
                1.0,
                pressure_scale * 101_325.0,
                pressure_scale * 101_325.0,
            ))
            .unwrap();
        assert!((result.temperature - reference.temperature).abs() <= 1e-6);
        assert!((result.candidate_moles - reference.candidate_moles).abs() <= 1e-9);
        assert!((result.total_enthalpy - reference.total_enthalpy).abs() <= 1e-7);
        for (actual, expected) in result.gas_moles.iter().zip(&reference.gas_moles) {
            assert!(
                (actual - expected).abs() <= 1e-9,
                "pressure_scale={pressure_scale:e}"
            );
        }
    }
}

#[test]
fn p10_1_ph_component_permutation_preserves_named_physical_state() {
    let validator = PurePhasePhValidator::default();
    let original_problem = synthetic_problem(1.0, 0.0);
    let permuted_problem = permuted_synthetic_problem();
    original_problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .expect("original family must remain structurally strict");
    permuted_problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
        .expect("permuted family must remain structurally strict");
    let original = validator.solve(&original_problem).unwrap();
    let permuted = validator.solve(&permuted_problem).unwrap();

    assert_eq!(permuted_problem.gas_species(), ["B(g)", "A(g)"]);
    assert!((permuted.temperature - original.temperature).abs() <= 1e-6);
    assert!((permuted.candidate_moles - original.candidate_moles).abs() <= 1e-9);
    assert!((permuted.total_enthalpy - original.total_enthalpy).abs() <= 1e-7);
    assert!((permuted.gas_moles[0] - original.gas_moles[1]).abs() <= 1e-9);
    assert!((permuted.gas_moles[1] - original.gas_moles[0]).abs() <= 1e-9);
}

#[test]
fn i1_i2_initially_present_candidate_uses_the_same_physical_coordinate() {
    // I1/I2: the pure candidate is already present, so the mathematical
    // extent is smaller while the recovered physical state is unchanged.
    let result = PurePhasePhValidator::default()
        .solve(&synthetic_problem_with_initial_candidate(1.0, 0.0, 0.10))
        .expect("an initially active pure phase remains a valid fixed-topology state");
    assert!((result.temperature - T_STAR).abs() <= 1e-6);
    assert!((result.extent - 0.15).abs() <= 1e-9);
    assert!((result.candidate_moles - XI_STAR).abs() <= 1e-9);
    assert!((result.gas_moles[0] - 1.5).abs() <= 1e-9);
    assert!((result.gas_moles[1] - 1.75).abs() <= 1e-9);
}

#[test]
fn p9_4_comparator_keeps_identity_and_missing_evidence_axes_explicit() {
    // P9.4: the full record supplies every comparable canonical axis.
    let problem = synthetic_problem(1.0, 0.0);
    let independent = PurePhasePhValidator::default().solve(&problem).unwrap();
    let evidence = PurePhasePhCanonicalEvidence::new(
        problem.case_identity(),
        independent.temperature,
        independent.gas_moles.clone(),
        independent.candidate_moles,
        independent.total_enthalpy,
        Some(independent.chemical_log_residual),
        Some(0.0),
    )
    .unwrap();
    let report = compare_pure_phase_ph_validation(
        &problem,
        &independent,
        &evidence,
        PurePhasePhCrossValidationTolerances::default(),
    )
    .unwrap();
    assert!(report.is_complete_match());

    let mismatched_case = synthetic_problem(1.0, 100.0);
    let mismatched_identity = PurePhasePhCanonicalEvidence::new(
        mismatched_case.case_identity(),
        independent.temperature,
        independent.gas_moles.clone(),
        independent.candidate_moles,
        independent.total_enthalpy,
        Some(0.0),
        Some(0.0),
    )
    .unwrap();
    let mismatch = compare_pure_phase_ph_validation(
        &problem,
        &independent,
        &mismatched_identity,
        PurePhasePhCrossValidationTolerances::default(),
    )
    .expect_err("canonical evidence from another P,H case must be rejected");
    assert!(matches!(
        mismatch,
        ReactionExtentError::InvalidProblem { .. }
    ));

    let foreign_result = PurePhasePhValidator::default()
        .solve(&synthetic_problem(1.0, 100.0))
        .expect("nearby synthetic target must remain independently solvable");
    let mismatch = compare_pure_phase_ph_validation(
        &problem,
        &foreign_result,
        &evidence,
        PurePhasePhCrossValidationTolerances::default(),
    )
    .expect_err("independent result from another target enthalpy must be rejected");
    assert!(matches!(
        mismatch,
        ReactionExtentError::InvalidProblem { .. }
    ));
    assert!(
        independent.max_abs_element_balance.is_some(),
        "strict synthetic I1/I2 result must publish independent conservation evidence"
    );

    // Matching identities alone do not turn unavailable canonical evidence
    // into a pass. This is the report contract used by partial real-data I4
    // comparisons when a particular accepted-state axis cannot be supplied.
    let partial = PurePhasePhCanonicalEvidence::new(
        problem.case_identity(),
        independent.temperature,
        independent.gas_moles.clone(),
        independent.candidate_moles,
        independent.total_enthalpy,
        None,
        None,
    )
    .unwrap();
    let partial_report = compare_pure_phase_ph_validation(
        &problem,
        &independent,
        &partial,
        PurePhasePhCrossValidationTolerances::default(),
    )
    .unwrap();
    assert_eq!(partial_report.identity_agreement, Some(true));
    assert_eq!(partial_report.thermodynamic_agreement, None);
    assert_eq!(partial_report.conservation_agreement, None);
    assert!(
        !partial_report.is_complete_match(),
        "missing evidence must never become an implicit successful comparison"
    );
}

#[test]
fn p9_4_comparator_localizes_each_physical_mismatch_axis() {
    // P9.4: preserve the identity, then perturb each physical evidence
    // axis independently. A failed cross-check must say where it disagreed.
    let problem = synthetic_problem(1.0, 0.0);
    let independent = PurePhasePhValidator::default().solve(&problem).unwrap();
    let evidence = PurePhasePhCanonicalEvidence::new(
        problem.case_identity(),
        independent.temperature + 1.0,
        vec![independent.gas_moles[0] + 0.1, independent.gas_moles[1]],
        independent.candidate_moles + 0.1,
        independent.total_enthalpy + 1.0,
        Some(1.0),
        Some(1.0),
    )
    .unwrap();
    let report = compare_pure_phase_ph_validation(
        &problem,
        &independent,
        &evidence,
        PurePhasePhCrossValidationTolerances::default(),
    )
    .unwrap();

    assert_eq!(report.identity_agreement, Some(true));
    assert_eq!(report.temperature_agreement, Some(false));
    assert_eq!(report.thermodynamic_agreement, Some(false));
    assert_eq!(report.composition_agreement, Some(false));
    assert_eq!(report.enthalpy_agreement, Some(false));
    assert_eq!(report.conservation_agreement, Some(false));
    assert!(!report.is_complete_match());
}

#[test]
fn p9_4_fixed_topology_monolithic_ph_matches_independent_nested_scalar_solution() {
    // P9.4 fixed-topology bridge, intentionally not I3 phase lifecycle.
    let problem = synthetic_problem_with_initial_candidate(1.0, 0.0, 0.1);
    let independent = PurePhasePhValidator::default().solve(&problem).unwrap();
    let canonical = canonical_fixed_topology_evidence(&problem);
    let report = compare_pure_phase_ph_validation(
        &problem,
        &independent,
        &canonical,
        PurePhasePhCrossValidationTolerances {
            max_abs_temperature_delta: 1e-4,
            max_abs_mole_delta: 1e-5,
            max_abs_enthalpy_delta: 1e-3,
            max_abs_chemical_log_residual: 1e-7,
            max_abs_element_balance: 1e-7,
            ..PurePhasePhCrossValidationTolerances::default()
        },
    )
    .unwrap();
    assert!(
        report.is_complete_match(),
        "P9 fixed-topology report: {report:?}"
    );
}

#[test]
fn i2_fixture_specific_enthalpy_ordering_is_monotone_near_constructed_state() {
    // I2: local fixture-specific H-target ordering around the constructed root.
    let validator = PurePhasePhValidator::default();
    // The thermodynamically consistent fixture has a compact outer enthalpy
    // range. Perturb inside its verified bracket rather than manufacturing an
    // unphysical extrapolation merely to assert an ordering.
    let lower = validator.solve(&synthetic_problem(1.0, -100.0)).unwrap();
    let center = validator.solve(&synthetic_problem(1.0, 0.0)).unwrap();
    let upper = validator.solve(&synthetic_problem(1.0, 100.0)).unwrap();
    assert!(lower.temperature < center.temperature);
    assert!(center.temperature < upper.temperature);
}

#[test]
fn i1_i2_failure_matrix_keeps_invalid_inputs_and_unbracketed_roots_typed() {
    // I1/I2: malformed domains and unverified scalar roots stay typed.
    assert!(matches!(
        PurePhasePhConditions::new(101_325.0, 101_325.0, 0.0, 900.0, 800.0),
        Err(ReactionExtentError::InvalidProblem { .. })
    ));
    assert!(matches!(
        PurePhasePhConditions::new(101_325.0, 101_325.0, 0.0, 0.0, 800.0),
        Err(ReactionExtentError::InvalidProblem { .. })
    ));
    assert!(matches!(
        PurePhasePhConditions::new(101_325.0, 101_325.0, f64::NAN, 650.0, 950.0),
        Err(ReactionExtentError::InvalidProblem { .. })
    ));

    let unbracketed = PurePhasePhValidator::default()
        .solve(&synthetic_problem(1.0, 1.0e8))
        .expect_err("a target outside the supplied enthalpy bracket must not be accepted");
    assert!(matches!(
        unbracketed,
        ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_ph_outer_temperature",
            ..
        }
    ));

    let no_inner_root = PurePhasePhProblem::new(
        vec!["A".into(), "B".into()],
        vec![2.0, 2.0],
        vec![-2.0, -1.0],
        0.0,
        1.0,
        "S",
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        function(|temperature| Ok(-MOLAR_GAS_CONSTANT * temperature * -100.0)),
        vec![function(|_| Ok(1.0)), function(|_| Ok(1.0))],
        function(|_| Ok(1.0)),
        PurePhasePhConditions::new(101_325.0, 101_325.0, 1.0, 650.0, 950.0).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        no_inner_root
            .solve_inner_extent_at_temperature(800.0, PurePhasePhSolverSettings::default()),
        Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_ph_inner_extent",
            ..
        })
    ));

    let exhausted = PurePhasePhValidator {
        settings: PurePhasePhSolverSettings {
            max_outer_iterations: 1,
            max_abs_enthalpy_residual: 1e-14,
            ..PurePhasePhSolverSettings::default()
        },
    }
    .solve(&synthetic_problem(1.0, 0.0))
    .expect_err("a one-step outer budget must not silently accept a midpoint");
    assert!(matches!(
        exhausted,
        ReactionExtentError::SolveError(SolveError::MaxIterations)
    ));

    let inner_budget_exhausted = synthetic_problem(1.0, 0.0)
        .solve_inner_extent_at_temperature(
            800.0,
            PurePhasePhSolverSettings {
                max_inner_iterations: 1,
                max_abs_log_residual: 1e-30,
                ..PurePhasePhSolverSettings::default()
            },
        )
        .expect_err("a one-step inner budget must not accept an unverified extent root");
    assert!(matches!(
        inner_budget_exhausted,
        ReactionExtentError::SolveError(SolveError::MaxIterations)
    ));

    let infeasible_extent = PurePhasePhProblem::new(
        vec!["A".into(), "B".into()],
        vec![2.0, 2.0],
        vec![1.0, 0.0],
        0.0,
        1.0,
        "S",
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        function(|_| Ok(0.0)),
        vec![function(|_| Ok(1.0)), function(|_| Ok(1.0))],
        function(|_| Ok(1.0)),
        PurePhasePhConditions::new(101_325.0, 101_325.0, 1.0, 650.0, 950.0).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        infeasible_extent
            .solve_inner_extent_at_temperature(800.0, PurePhasePhSolverSettings::default()),
        Err(ReactionExtentError::ValidationNotApplicable {
            path: "pure_phase_ph_inner_extent",
            ..
        })
    ));
}

#[test]
fn p10_1_ph_settings_and_comparison_tolerances_reject_invalid_values() {
    let problem = synthetic_problem(1.0, 0.0);
    let independent = PurePhasePhValidator::default()
        .solve(&problem)
        .expect("baseline P,H scalar problem must solve");
    let canonical = PurePhasePhCanonicalEvidence::new(
        problem.case_identity(),
        independent.temperature,
        independent.gas_moles.clone(),
        independent.candidate_moles,
        independent.total_enthalpy,
        Some(independent.chemical_log_residual),
        independent.max_abs_element_balance,
    )
    .expect("baseline P,H canonical evidence must validate");

    for settings in [
        PurePhasePhSolverSettings {
            max_inner_iterations: 0,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            max_outer_iterations: 0,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            interior_branch_scan_subdivisions: 1,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            max_abs_log_residual: f64::NAN,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            max_abs_enthalpy_residual: -1.0,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            max_abs_enthalpy_residual: f64::INFINITY,
            ..PurePhasePhSolverSettings::default()
        },
        PurePhasePhSolverSettings {
            feasibility_margin: 0.25,
            ..PurePhasePhSolverSettings::default()
        },
    ] {
        let error = PurePhasePhValidator { settings }
            .solve(&problem)
            .expect_err("invalid P,H scalar settings must be rejected before solving");
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }

    for tolerances in [
        PurePhasePhCrossValidationTolerances {
            max_abs_mole_delta: 0.0,
            ..PurePhasePhCrossValidationTolerances::default()
        },
        PurePhasePhCrossValidationTolerances {
            max_relative_temperature_delta: -1.0,
            ..PurePhasePhCrossValidationTolerances::default()
        },
        PurePhasePhCrossValidationTolerances {
            max_relative_mole_delta: f64::NAN,
            ..PurePhasePhCrossValidationTolerances::default()
        },
        PurePhasePhCrossValidationTolerances {
            max_relative_enthalpy_delta: f64::INFINITY,
            ..PurePhasePhCrossValidationTolerances::default()
        },
    ] {
        let error =
            compare_pure_phase_ph_validation(&problem, &independent, &canonical, tolerances)
                .expect_err("invalid P,H comparison tolerance must be rejected before comparison");
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }

    for structural in [
        PurePhaseBoundaryStructuralTolerances {
            rank_absolute_tolerance: 0.0,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
        PurePhaseBoundaryStructuralTolerances {
            max_abs_element_balance: -1.0,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
        PurePhaseBoundaryStructuralTolerances {
            rank_relative_tolerance: f64::NAN,
            ..PurePhaseBoundaryStructuralTolerances::default()
        },
    ] {
        let error = problem
            .reaction_space(structural)
            .expect_err("invalid P,H structural tolerance must be rejected");
        assert!(matches!(error, ReactionExtentError::InvalidProblem { .. }));
    }
}

#[test]
fn i1_i2_retargeted_enthalpy_problem_is_immutable_and_revalidated() {
    // I1/I2: target-range validation must not mutate the reference problem
    // retained for a preceding accepted point.
    let original = synthetic_problem(1.0, 0.0);
    let original_target = original.conditions().target_enthalpy();
    let retargeted = original
        .with_target_enthalpy(original_target + 100.0)
        .unwrap();
    assert_eq!(original.conditions().target_enthalpy(), original_target);
    assert_eq!(
        retargeted.conditions().target_enthalpy(),
        original_target + 100.0
    );
    assert!(
        PurePhasePhValidator::default().solve(&retargeted).is_ok(),
        "a neighbouring fixture target must retain a valid scalar reference"
    );
    assert!(matches!(
        original.with_target_enthalpy(f64::NAN),
        Err(ReactionExtentError::InvalidProblem { .. })
    ));
}

#[test]
fn i1_i2_non_finite_thermochemistry_is_not_converted_to_nan() {
    // I1/I2: fallible thermo capabilities cannot leak NaN into scalar solving.
    let problem = PurePhasePhProblem::new(
        vec!["A".into(), "B".into()],
        vec![2.0, 2.0],
        vec![-2.0, -1.0],
        0.0,
        1.0,
        "S",
        vec![function(|_| Ok(f64::NAN)), function(|_| Ok(0.0))],
        function(|_| Ok(0.0)),
        vec![function(|_| Ok(1.0)), function(|_| Ok(1.0))],
        function(|_| Ok(1.0)),
        PurePhasePhConditions::new(101_325.0, 101_325.0, 1.0, 650.0, 950.0).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        PurePhasePhValidator::default().solve(&problem),
        Err(ReactionExtentError::InvalidDG0 {
            species_index: 0,
            ..
        })
    ));

    let invalid_enthalpy = PurePhasePhProblem::new(
        vec!["A".into(), "B".into()],
        vec![2.0, 2.0],
        vec![-2.0, -1.0],
        0.0,
        1.0,
        "S",
        vec![function(|_| Ok(0.0)), function(|_| Ok(0.0))],
        function(|temperature| Ok(-MOLAR_GAS_CONSTANT * temperature * 3.0)),
        vec![function(|_| Ok(1.0)), function(|_| Ok(1.0))],
        function(|_| Ok(f64::INFINITY)),
        PurePhasePhConditions::new(101_325.0, 101_325.0, 1.0, 650.0, 950.0).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        PurePhasePhValidator::default().solve(&invalid_enthalpy),
        Err(ReactionExtentError::InvalidCandidate {
            field: "pure_phase_ph_enthalpy",
            ..
        })
    ));
}
