use super::*;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_calculator::{
    EquilibriumCalculator, EquilibriumCalculatorOutcome, EquilibriumCalculatorPreset,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::{
    RstPreparedProblem, RstSymbolicFrontend, RustedSciTheSolver,
    prepare_rst_symbolic_ph_problem_with_timing,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_argonne_stanjan_chon::STANJAN_LOCAL_GAS_SPECIES;
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use RustedSciThe::numerical::Nonlinear_systems::prelude::{
    LambdifyExecutionPolicy, NonlinearProblem, PreparationTelemetryMode,
    PreparedSymbolicNonlinearProblem, SymbolicLambdifyFrontend, SymbolicProblemOptions,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::fmt::Write as _;

const TEMPERATURE_K: f64 = 2_500.0;
const PRESSURE_PA: f64 = 101_325.0;
const REFERENCE_PRESSURE_PA: f64 = 100_000.0;

fn seed_moles() -> Vec<f64> {
    STANJAN_LOCAL_GAS_SPECIES
        .iter()
        .map(|species| match *species {
            "CH4" => 4.0,
            "CO2" => 2.0,
            "O2" => 8.0,
            "N2" => 37.6,
            _ => 0.0,
        })
        .collect()
}

fn diagnostic_inputs() -> (ResolvedPhaseEnthalpyRequest<'static>, PreparedPhFormulation) {
    let phase = PhaseId::new(Some("gas".to_string()));
    let initial = seed_moles();
    let pt_outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", STANJAN_LOCAL_GAS_SPECIES)
        .expect("STANJAN gas phase")
        .initial_phase_moles(
            STANJAN_LOCAL_GAS_SPECIES
                .iter()
                .zip(initial.iter().copied())
                .filter(|(_, moles)| *moles > 0.0)
                .map(|(species, moles)| (PhaseComponentId::new(phase.clone(), *species), moles)),
        )
        .prefer_libraries(["NASA_gas"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(REFERENCE_PRESSURE_PA)
        .at_temperature(TEMPERATURE_K)
        .production_cascade()
        .preset(EquilibriumCalculatorPreset::Silent)
        .solve()
        .expect("STANJAN reference P,T solve");
    let EquilibriumCalculatorOutcome::PtPoint(pt_point) = pt_outcome else {
        panic!("expected STANJAN P,T point")
    };
    let resolved = pt_point.resolved().clone();
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("resolved STANJAN thermochemistry");
    let target_enthalpy = thermochemistry
        .enthalpy_model()
        .evaluate_total(pt_point.solution().component_moles(), TEMPERATURE_K)
        .expect("reference total enthalpy");

    let layout =
        MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec()).expect("STANJAN layout");
    let initial_composition =
        MultiphaseInitialComposition::from_dense(&layout, initial).expect("STANJAN molecular seed");
    let temperature_bounds =
        TemperatureBounds::new(2_100.0, 2_900.0).expect("STANJAN P,H temperature bounds");
    let constraint = EquilibriumConstraint::ph(
        PRESSURE_PA,
        REFERENCE_PRESSURE_PA,
        target_enthalpy,
        TEMPERATURE_K,
    )
    .expect("STANJAN P,H constraint");
    let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
        &resolved,
        initial_composition,
        constraint,
        temperature_bounds,
        thermochemistry,
    )
    .expect("resolved STANJAN P,H request")
    .with_ph_solve_mode(PhSolveMode::Monolithic);

    let (target, initial_temperature) =
        validated_ph_parameters(request.constraint).expect("validated P,H target");
    let thermo = request
        .thermochemistry
        .clone()
        .expect("resolved thermochemistry on request");
    let initial_enthalpies = thermo
        .evaluate_enthalpy(initial_temperature)
        .expect("initial enthalpy values");
    let enthalpy_scale = EnthalpyScale::from_magnitudes(
        target,
        request.initial_composition.moles(),
        &initial_enthalpies,
    )
    .expect("P,H enthalpy scale");
    let conditions = request
        .constraint
        .conditions_at(initial_temperature)
        .expect("seed P,T conditions");
    let prepared = build_phase_equilibrium_problem_with_timing(
        request
            .build_pt_request(conditions)
            .expect("fixed P,T build request"),
        EquilibriumTimingMode::Disabled,
    )
    .expect("fixed P,T preparation")
    .into_prepared_fixed_active_problem()
    .expect("fixed active P,T problem");
    let formulation = PreparedPhFormulation::new(
        prepared.prepared,
        thermo,
        temperature_bounds,
        target,
        enthalpy_scale,
    )
    .expect("native monolithic P,H formulation");

    (request, formulation)
}

fn prepare_frontend(
    formulation: &PreparedPhFormulation,
    frontend: RstSymbolicFrontend,
) -> RstPreparedProblem {
    prepare_rst_symbolic_ph_problem_with_timing(
        formulation,
        EquilibriumTimingMode::Disabled,
        frontend,
    )
    .unwrap_or_else(|error| panic!("{frontend:?} preparation failed: {error}"))
    .0
}

fn max_abs_vector_delta(left: &[f64], right: &[f64]) -> (f64, usize) {
    assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .enumerate()
        .map(|(index, (left, right))| (index, (left - right).abs()))
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map_or((0.0, 0), |(index, delta)| (delta, index))
}

fn max_abs_matrix_delta(
    left: &nalgebra::DMatrix<f64>,
    right: &nalgebra::DMatrix<f64>,
) -> (f64, (usize, usize)) {
    assert_eq!(left.shape(), right.shape());
    let rows = left.nrows().max(1);
    left.iter()
        .zip(right.iter())
        .enumerate()
        .map(|(flat, (left, right))| {
            let row = flat % rows;
            let col = flat / rows;
            ((row, col), (left - right).abs())
        })
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map_or((0.0, (0, 0)), |(index, delta)| (delta, index))
}

fn epsilon_units(left: f64, right: f64) -> f64 {
    (left - right).abs() / (f64::EPSILON * left.abs().max(right.abs()).max(1.0))
}

fn expression_children(expression: &Expr) -> Vec<&Expr> {
    use Expr::*;
    match expression {
        Var(_) | Const(_) => Vec::new(),
        Add(left, right)
        | Sub(left, right)
        | Mul(left, right)
        | Div(left, right)
        | Pow(left, right) => vec![left, right],
        Exp(value) | Ln(value) | sin(value) | cos(value) | tg(value) | ctg(value)
        | arcsin(value) | arccos(value) | arctg(value) | arcctg(value) => vec![value],
    }
}

fn evaluate_expression_frontend(
    expression: &Expr,
    frontend: SymbolicLambdifyFrontend,
    variables: &[String],
    unknowns: &[f64],
    parameters: &[String],
    parameter_values: &[f64],
) -> f64 {
    let mut equations = Vec::with_capacity(variables.len());
    equations.push(expression.clone());
    equations.resize_with(variables.len(), || Expr::Const(0.0));
    let mut options = SymbolicProblemOptions::new()
        .with_variables(variables.to_vec())
        .with_lambdify_execution_policy(LambdifyExecutionPolicy::Sequential)
        .with_lambdify_frontend(frontend)
        .with_preparation_telemetry(PreparationTelemetryMode::Disabled);
    if !parameters.is_empty() {
        options = options
            .with_equation_parameters(parameters.to_vec())
            .with_equation_parameter_values(nalgebra::DVector::from_column_slice(parameter_values));
    }
    let prepared = PreparedSymbolicNonlinearProblem::from_expressions(equations, options)
        .unwrap_or_else(|error| panic!("{frontend:?} subtree preparation failed: {error}"));
    let bound = prepared
        .bind_initial()
        .unwrap_or_else(|error| panic!("{frontend:?} subtree binding failed: {error}"));
    let result = if parameters.is_empty() {
        bound
            .prepared()
            .as_problem()
            .residual(&nalgebra::DVector::from_column_slice(unknowns))
            .unwrap_or_else(|error| panic!("{frontend:?} subtree evaluation failed: {error}"))
    } else {
        bound
            .residual_with_parameter_values(
                &nalgebra::DVector::from_column_slice(unknowns),
                &nalgebra::DVector::from_column_slice(parameter_values),
            )
            .unwrap_or_else(|error| panic!("{frontend:?} subtree evaluation failed: {error}"))
    };
    result[0]
}

fn shortened_expression(expression: &Expr) -> String {
    let text = expression.to_string();
    const LIMIT: usize = 320;
    if text.chars().count() <= LIMIT {
        text
    } else {
        format!(
            "{}... <{} chars>",
            text.chars().take(LIMIT).collect::<String>(),
            text.chars().count()
        )
    }
}

fn report_subexpression_descent(
    report: &mut String,
    root: &Expr,
    variables: &[String],
    unknowns: &[f64],
    parameters: &[String],
    parameter_values: &[f64],
) {
    writeln!(
        report,
        "\nsubexpression descent (same coordinates and bound parameters):"
    )
    .unwrap();
    let mut current = root;
    for depth in 0..64 {
        let expr_value = evaluate_expression_frontend(
            current,
            SymbolicLambdifyFrontend::ExprLegacy,
            variables,
            unknowns,
            parameters,
            parameter_values,
        );
        let atom_value = evaluate_expression_frontend(
            current,
            SymbolicLambdifyFrontend::AtomViewNative,
            variables,
            unknowns,
            parameters,
            parameter_values,
        );
        let delta = (expr_value - atom_value).abs();
        writeln!(
            report,
            "  depth={depth} node ExprLegacy={expr_value:.17e} AtomView={atom_value:.17e} abs_delta={delta:.17e} expr={}",
            shortened_expression(current),
        )
        .unwrap();

        let children = expression_children(current);
        if children.is_empty() {
            writeln!(report, "    stop=leaf").unwrap();
            break;
        }
        let mut scored = Vec::with_capacity(children.len());
        for (index, child) in children.iter().enumerate() {
            let expr_child = evaluate_expression_frontend(
                child,
                SymbolicLambdifyFrontend::ExprLegacy,
                variables,
                unknowns,
                parameters,
                parameter_values,
            );
            let atom_child = evaluate_expression_frontend(
                child,
                SymbolicLambdifyFrontend::AtomViewNative,
                variables,
                unknowns,
                parameters,
                parameter_values,
            );
            let child_delta = (expr_child - atom_child).abs();
            let score = child_delta / expr_child.abs().max(atom_child.abs()).max(1.0);
            writeln!(
                report,
                "    child[{index}] ExprLegacy={expr_child:.17e} AtomView={atom_child:.17e} abs_delta={child_delta:.17e} normalized_delta={score:.17e} expr={}",
                shortened_expression(child),
            )
            .unwrap();
            scored.push((index, score));
        }
        let Some((child_index, score)) = scored
            .into_iter()
            .max_by(|left, right| left.1.total_cmp(&right.1))
        else {
            break;
        };
        if score <= 1.0e-14 {
            if delta > 1.0e-12 * expr_value.abs().max(atom_value.abs()).max(1.0) {
                writeln!(
                    report,
                    "    stop=parent_discrepancy_is_amplified_from_sub_roundoff_child_delta; parent_expr_exact={current}"
                )
                .unwrap();
            } else {
                writeln!(report, "    stop=children_agree_to_roundoff").unwrap();
            }
            break;
        }
        writeln!(report, "    descend=child[{child_index}]").unwrap();
        current = children[child_index];
    }
}

fn report_tiny_constant_lambdify_probe() {
    const TEMPERATURE: f64 = 2_500.0;
    const COEFFICIENT_TEXT: &str = "-0.000000000000000509024475";
    let coefficient = Expr::parse_expression(COEFFICIENT_TEXT);
    let polynomial_term = Expr::parse_expression("-0.000000000000000509024475 * T ^ 4");
    let variables = vec!["T".to_string()];
    let atom_coefficient = evaluate_expression_frontend(
        &coefficient,
        SymbolicLambdifyFrontend::AtomViewNative,
        &variables,
        &[TEMPERATURE],
        &[],
        &[],
    );
    let atom_term = evaluate_expression_frontend(
        &polynomial_term,
        SymbolicLambdifyFrontend::AtomViewNative,
        &variables,
        &[TEMPERATURE],
        &[],
        &[],
    );
    let variable_names = ["T"];
    let direct_expr_coefficient = coefficient.lambdify1(&variable_names)(&[TEMPERATURE]);
    let direct_expr_term = polynomial_term.lambdify1(&variable_names)(&[TEMPERATURE]);
    let prepared_expr_coefficient = evaluate_expression_frontend(
        &coefficient,
        SymbolicLambdifyFrontend::ExprLegacy,
        &variables,
        &[TEMPERATURE],
        &[],
        &[],
    );
    let prepared_expr_term = evaluate_expression_frontend(
        &polynomial_term,
        SymbolicLambdifyFrontend::ExprLegacy,
        &variables,
        &[TEMPERATURE],
        &[],
        &[],
    );
    let within_f64_roundoff = |left: f64, right: f64| {
        (left - right).abs() <= 8.0 * f64::EPSILON * left.abs().max(right.abs())
    };

    let mut report = String::new();
    writeln!(report, "Minimal ExprLegacy versus AtomView Lambdify probe").unwrap();
    writeln!(report, "T={TEMPERATURE:.17e}").unwrap();
    writeln!(report, "coefficient_source={COEFFICIENT_TEXT}").unwrap();
    writeln!(report, "coefficient parsed Expr={coefficient}").unwrap();
    writeln!(
        report,
        "coefficient direct_Expr_lambdify1={direct_expr_coefficient:.17e} prepared_ExprLegacy={prepared_expr_coefficient:.17e} AtomView={atom_coefficient:.17e}"
    )
    .unwrap();
    writeln!(report, "term_source=-0.000000000000000509024475 * T ^ 4").unwrap();
    writeln!(report, "term parsed Expr={polynomial_term}").unwrap();
    writeln!(
        report,
        "term direct_Expr_lambdify1={direct_expr_term:.17e} prepared_ExprLegacy={prepared_expr_term:.17e} AtomView={atom_term:.17e}"
    )
    .unwrap();
    writeln!(report, "T^4={:.17e}", TEMPERATURE.powi(4)).unwrap();
    writeln!(report, "coefficient_as_f64={:.17e}", -5.09024475e-16_f64).unwrap();

    for coefficient_text in [
        "0.00000000000509024475",
        "0.000000000000509024475",
        "0.0000000000000509024475",
        "0.00000000000000509024475",
        "0.000000000000000509024475",
        "0.0000000000000000509024475",
        "0.00000000000000000509024475",
        "-0.00000000000509024475",
        "-0.000000000000509024475",
        "-0.0000000000000509024475",
        "-0.00000000000000509024475",
        COEFFICIENT_TEXT,
        "-0.0000000000000000509024475",
        "-0.00000000000000000509024475",
    ] {
        let parsed = Expr::parse_expression(coefficient_text);
        let expr_value = parsed.lambdify1(&variable_names)(&[TEMPERATURE]);
        let atom_value = evaluate_expression_frontend(
            &parsed,
            SymbolicLambdifyFrontend::AtomViewNative,
            &variables,
            &[TEMPERATURE],
            &[],
            &[],
        );
        assert!(
            within_f64_roundoff(atom_value, expr_value),
            "AtomView constant mismatch for {coefficient_text}: ExprLegacy={expr_value:.17e}, AtomView={atom_value:.17e}"
        );
        writeln!(
            report,
            "constant_sweep source={coefficient_text} ExprLegacy={expr_value:.17e} AtomView={atom_value:.17e} abs_delta={:.17e}",
            (expr_value - atom_value).abs(),
        )
        .unwrap();
    }

    assert_ne!(direct_expr_coefficient, 0.0);
    assert_eq!(prepared_expr_coefficient, direct_expr_coefficient);
    assert!(
        within_f64_roundoff(atom_coefficient, direct_expr_coefficient),
        "AtomView lost the parsed nonzero coefficient: ExprLegacy={direct_expr_coefficient:.17e}, AtomView={atom_coefficient:.17e}"
    );
    assert_eq!(prepared_expr_term, direct_expr_term);
    assert!(
        within_f64_roundoff(atom_term, direct_expr_term),
        "AtomView changed the amplified polynomial term: ExprLegacy={direct_expr_term:.17e}, AtomView={atom_term:.17e}"
    );
    let path = "target/stanjan15_atom_tiny_constant_probe.txt";
    std::fs::write(path, report).expect("write tiny-constant Lambdify probe");
    println!(
        "tiny_constant_probe ExprLegacy_coef={direct_expr_coefficient:.17e} AtomView_coef={atom_coefficient:.17e} ExprLegacy_term={direct_expr_term:.17e} AtomView_term={atom_term:.17e} report={path}"
    );
}

fn central_difference(
    problem: &RstPreparedProblem,
    unknowns: &[f64],
    row: usize,
    column: usize,
) -> f64 {
    let step = f64::EPSILON.cbrt() * unknowns[column].abs().max(1.0);
    let mut plus = unknowns.to_vec();
    let mut minus = unknowns.to_vec();
    plus[column] += step;
    minus[column] -= step;
    let plus_value = problem
        .residual_for_test(&plus)
        .expect("positive finite-difference residual")[row];
    let minus_value = problem
        .residual_for_test(&minus)
        .expect("negative finite-difference residual")[row];
    (plus_value - minus_value) / (2.0 * step)
}

fn report_symbolic_witness(
    label: &str,
    expr: &RstPreparedProblem,
    atom: &RstPreparedProblem,
    unknowns: &[f64],
    residual_row: usize,
    jacobian_index: (usize, usize),
    expr_residual: &[f64],
    atom_residual: &[f64],
    native_residual: &[f64],
    expr_jacobian: &nalgebra::DMatrix<f64>,
    atom_jacobian: &nalgebra::DMatrix<f64>,
    native_jacobian: &nalgebra::DMatrix<f64>,
) {
    let (expr_equations, variables, parameter_values, parameters) =
        expr.symbolic_payload_for_test();
    let (atom_equations, atom_variables, atom_parameter_values, atom_parameters) =
        atom.symbolic_payload_for_test();
    let expr_jacobian_expressions = expr.symbolic_jacobian_for_test();
    let atom_jacobian_expressions = atom.symbolic_jacobian_for_test();
    assert_eq!(variables, atom_variables);
    assert_eq!(parameters, atom_parameters);
    assert_eq!(parameter_values, atom_parameter_values);

    let mut report = String::new();
    writeln!(
        report,
        "STANJAN-15 P,H AtomView/ExprLegacy witness: {label}"
    )
    .unwrap();
    writeln!(report, "coordinates=[ln(moles), transformed_temperature]").unwrap();
    writeln!(
        report,
        "equation_strings_identical={}",
        expr_equations == atom_equations
    )
    .unwrap();
    writeln!(report, "x coordinates:").unwrap();
    for (index, (name, value)) in variables.iter().zip(unknowns).enumerate() {
        writeln!(report, "  x[{index}] {name} = {value:.17e}").unwrap();
    }
    for (index, (name, value)) in parameters.iter().zip(&parameter_values).enumerate() {
        writeln!(report, "  parameter[{index}] {name} = {value:.17e}").unwrap();
    }
    writeln!(report, "\nresidual equations and evaluations:").unwrap();
    for row in 0..expr_residual.len() {
        let delta = (expr_residual[row] - atom_residual[row]).abs();
        writeln!(
            report,
            "  row={row} ExprLegacy={:.17e} AtomView={:.17e} native={:.17e} abs_delta={delta:.17e} epsilon_units={:.6e}",
            expr_residual[row],
            atom_residual[row],
            native_residual[row],
            epsilon_units(expr_residual[row], atom_residual[row]),
        )
        .unwrap();
    }
    writeln!(report, "\nExprLegacy symbolic residual[{residual_row}]:").unwrap();
    writeln!(report, "{}", expr_equations[residual_row]).unwrap();
    if expr_equations[residual_row] != atom_equations[residual_row] {
        writeln!(report, "\nAtomView source residual[{residual_row}]:").unwrap();
        writeln!(report, "{}", atom_equations[residual_row]).unwrap();
    }
    let symbolic_equations = expr.symbolic_equations_for_test();
    report_subexpression_descent(
        &mut report,
        &symbolic_equations[residual_row],
        &variables,
        unknowns,
        &parameters,
        &parameter_values,
    );
    let row = jacobian_index.0;
    let column = jacobian_index.1;
    let finite_expr = central_difference(expr, unknowns, row, column);
    let finite_atom = central_difference(atom, unknowns, row, column);
    writeln!(
        report,
        "\nJacobian entries with |ExprLegacy - AtomView| > 1e-10 (numeric values):"
    )
    .unwrap();
    for i in 0..expr_jacobian.nrows() {
        for j in 0..expr_jacobian.ncols() {
            let expr_value = expr_jacobian[(i, j)];
            let atom_value = atom_jacobian[(i, j)];
            let delta = (expr_value - atom_value).abs();
            if delta <= 1.0e-10 {
                continue;
            }
            writeln!(
                report,
                "  J[{i},{j}] d/d{} ExprLegacy={expr_value:.17e} AtomView={atom_value:.17e} native={:.17e} abs_delta={delta:.17e} epsilon_units={:.6e}",
                variables[j],
                native_jacobian[(i, j)],
                epsilon_units(expr_value, atom_value),
            )
            .unwrap();
        }
    }
    writeln!(
        report,
        "\nselected Jacobian witness J[{row},{column}] = d residual[{row}] / d {}:",
        variables[column]
    )
    .unwrap();
    writeln!(report, "  residual_expr={}", expr_equations[row]).unwrap();
    if let Some(expression) = expr_jacobian_expressions
        .get(row)
        .and_then(|line| line.get(column))
    {
        writeln!(report, "  derivative_expr_ExprLegacy={expression}").unwrap();
    }
    if let Some(expression) = atom_jacobian_expressions
        .get(row)
        .and_then(|line| line.get(column))
    {
        writeln!(report, "  derivative_expr_AtomView={expression}").unwrap();
    } else {
        writeln!(
            report,
            "  derivative_expr_AtomView=<not materialized by AtomView; its native Jacobian value is reported below>"
        )
        .unwrap();
    }
    writeln!(
        report,
        "  evaluated ExprLegacy={:.17e} AtomView={:.17e} native={:.17e} abs_delta={:.17e} epsilon_units={:.6e}",
        expr_jacobian[(row, column)],
        atom_jacobian[(row, column)],
        native_jacobian[(row, column)],
        (expr_jacobian[(row, column)] - atom_jacobian[(row, column)]).abs(),
        epsilon_units(expr_jacobian[(row, column)], atom_jacobian[(row, column)]),
    )
    .unwrap();
    writeln!(
        report,
        "  central_difference ExprLegacy={finite_expr:.17e} AtomView={finite_atom:.17e}"
    )
    .unwrap();

    let path = format!("target/stanjan15_atom_view_witness_{label}.txt");
    std::fs::write(&path, report).expect("write detailed AtomView witness report");
    println!(
        "symbolic_witness stage={label} max_residual_delta={:.9e} row={residual_row} max_jacobian_delta={:.9e} index=({row},{column}) report={path}",
        (expr_residual[residual_row] - atom_residual[residual_row]).abs(),
        (expr_jacobian[(row, column)] - atom_jacobian[(row, column)]).abs(),
    );
    if label == "initial_seed" {
        println!("ExprLegacy/AtomView shared source expression residual[{residual_row}]:");
        println!("{}", expr_equations[residual_row]);
        println!("ExprLegacy symbolic expression for Jacobian[{row},{column}]:");
        if let Some(expression) = expr_jacobian_expressions
            .get(row)
            .and_then(|line| line.get(column))
        {
            println!("{expression}");
        }
        println!(
            "AtomView has no retained Expr Jacobian; compare the AtomView callback value and central-difference witness in {path}"
        );
    }
}

fn compare_at(
    label: &str,
    formulation: &PreparedPhFormulation,
    expr: &RstPreparedProblem,
    atom: &RstPreparedProblem,
    unknowns: &[f64],
) {
    let native_residual = formulation
        .residual(unknowns)
        .unwrap_or_else(|error| panic!("native {label} residual failed: {error}"));
    let native_jacobian = formulation
        .jacobian(unknowns)
        .unwrap_or_else(|error| panic!("native {label} Jacobian failed: {error}"));
    let expr_residual = expr
        .residual_for_test(unknowns)
        .unwrap_or_else(|error| panic!("ExprLegacy {label} residual failed: {error}"));
    let atom_residual = atom
        .residual_for_test(unknowns)
        .unwrap_or_else(|error| panic!("AtomView {label} residual failed: {error}"));
    let expr_jacobian = expr
        .jacobian_for_test(unknowns)
        .unwrap_or_else(|error| panic!("ExprLegacy {label} Jacobian failed: {error}"));
    let atom_jacobian = atom
        .jacobian_for_test(unknowns)
        .unwrap_or_else(|error| panic!("AtomView {label} Jacobian failed: {error}"));
    let (residual_delta, residual_row) =
        max_abs_vector_delta(expr_residual.as_slice(), atom_residual.as_slice());
    let (jacobian_delta, jacobian_index) = max_abs_matrix_delta(&expr_jacobian, &atom_jacobian);
    assert!(
        residual_delta <= 1.0e-12,
        "ExprLegacy/AtomView residual parity failed at {label}: max delta {residual_delta:.17e} at row {residual_row}"
    );
    assert!(
        jacobian_delta <= 1.0e-12,
        "ExprLegacy/AtomView Jacobian parity failed at {label}: max delta {jacobian_delta:.17e} at {jacobian_index:?}"
    );
    let (atom_native_residual_delta, atom_native_residual_row) =
        max_abs_vector_delta(atom_residual.as_slice(), &native_residual);
    let (atom_native_jacobian_delta, atom_native_jacobian_index) =
        max_abs_matrix_delta(&atom_jacobian, &native_jacobian);
    println!(
        "linearization stage={label} dimensions={} expr_residual_l2={:.9e} atom_residual_l2={:.9e} native_residual_l2={:.9e} expr_atom_residual_max={:.9e}@row{} atom_native_residual_max={:.9e}@row{} expr_atom_jacobian_max={:.9e}@({},{}) atom_native_jacobian_max={:.9e}@({},{})",
        unknowns.len(),
        expr_residual.norm(),
        atom_residual.norm(),
        nalgebra::DVector::from_column_slice(&native_residual).norm(),
        residual_delta,
        residual_row,
        atom_native_residual_delta,
        atom_native_residual_row,
        jacobian_delta,
        jacobian_index.0,
        jacobian_index.1,
        atom_native_jacobian_delta,
        atom_native_jacobian_index.0,
        atom_native_jacobian_index.1,
    );
    println!(
        "  worst residual row {}: expr={:.17e} atom={:.17e} native={:.17e}; worst Jacobian ({},{}): expr={:.17e} atom={:.17e} native={:.17e}",
        residual_row,
        expr_residual[residual_row],
        atom_residual[residual_row],
        native_residual[residual_row],
        jacobian_index.0,
        jacobian_index.1,
        expr_jacobian[(jacobian_index.0, jacobian_index.1)],
        atom_jacobian[(jacobian_index.0, jacobian_index.1)],
        native_jacobian[(jacobian_index.0, jacobian_index.1)],
    );
    if label == "initial_seed" || label == "AtomView" {
        report_symbolic_witness(
            label,
            expr,
            atom,
            unknowns,
            residual_row,
            jacobian_index,
            expr_residual.as_slice(),
            atom_residual.as_slice(),
            &native_residual,
            &expr_jacobian,
            &atom_jacobian,
            &native_jacobian,
        );
    }
}

fn solve_with_capture(
    request: &ResolvedPhaseEnthalpyRequest<'_>,
    frontend: RstSymbolicFrontend,
) -> Vec<(Vec<f64>, usize)> {
    let options = EquilibriumSolveOptions::new()
        .with_solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
            RustedSciTheSolver::MinpackLevenbergMarquardt,
        )))
        .expect("single-Minpack policy")
        .with_rst_symbolic_frontend(frontend);
    let request = request
        .clone()
        .with_solve_options(options)
        .with_ph_solve_mode(PhSolveMode::Monolithic);
    RustedSciTheSolver::begin_history_capture_for_test();
    let outcome = solve_resolved_ph(request);
    let captures = RustedSciTheSolver::take_history_capture_for_test();
    println!(
        "facade frontend={frontend:?} accepted={} result={:?} captured_attempts={}",
        outcome.is_ok(),
        outcome.as_ref().err(),
        captures.len()
    );
    let mut finals = Vec::with_capacity(captures.len());
    for (attempt, result) in captures.iter().enumerate() {
        println!(
            "attempt frontend={frontend:?} index={attempt} termination={:?} iterations={} residual_l2={:.9e} trace_len={}",
            result.termination,
            result.iterations,
            result.residual_norm,
            result.history.len(),
        );
        for record in result
            .history
            .iter()
            .take(4)
            .chain(result.history.iter().rev().take(4).rev())
        {
            println!(
                "  iteration={} accepted={} residual_l2={:.9e} step_norm={:.9e}",
                record.iteration, record.accepted, record.residual_norm, record.step_norm,
            );
        }
        finals.push((result.x.as_slice().to_vec(), result.history.len()));
    }
    finals
}

fn solve_public_ph_facade(
    target_enthalpy_j: f64,
    frontend: RstSymbolicFrontend,
) -> (f64, Vec<f64>, usize) {
    let phase = PhaseId::new(Some("gas".to_string()));
    let initial = seed_moles();
    let outcome = EquilibriumCalculator::builder()
        .ideal_gas_phase("gas", STANJAN_LOCAL_GAS_SPECIES)
        .expect("STANJAN facade phase")
        .initial_phase_moles(
            STANJAN_LOCAL_GAS_SPECIES
                .iter()
                .zip(initial.iter().copied())
                .filter(|(_, moles)| *moles > 0.0)
                .map(|(species, moles)| (PhaseComponentId::new(phase.clone(), *species), moles)),
        )
        .prefer_libraries(["NASA_gas"])
        .offline_only()
        .pressure_pa(PRESSURE_PA)
        .reference_pressure_pa(REFERENCE_PRESSURE_PA)
        .at_total_enthalpy(
            crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TotalEnthalpyJoules::new(
                target_enthalpy_j,
            )
            .expect("finite target enthalpy"),
            TEMPERATURE_K,
            TemperatureBounds::new(2_100.0, 2_900.0).expect("temperature bounds"),
        )
        .ph_solve_mode(PhSolveMode::Monolithic)
        .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
            RustedSciTheSolver::MinpackLevenbergMarquardt,
        )))
        .expect("single-Minpack facade policy")
        .rst_symbolic_frontend(frontend)
        .preset(EquilibriumCalculatorPreset::Silent)
        .solve()
        .unwrap_or_else(|error| panic!("public P,H facade {frontend:?} failed: {error}"));
    let EquilibriumCalculatorOutcome::PhPoint(point) = outcome else {
        panic!("expected public P,H point for {frontend:?}")
    };
    (
        point.solution().temperature(),
        point.solution().equilibrium().component_moles().to_vec(),
        point.solution().equilibrium().phases().len(),
    )
}

#[test]
fn stanjan15_tiny_constant_and_public_ph_facade_are_frontend_equivalent() {
    report_tiny_constant_lambdify_probe();
    let (request, formulation) = diagnostic_inputs();
    let expr = prepare_frontend(&formulation, RstSymbolicFrontend::ExprLegacy);
    let atom = prepare_frontend(&formulation, RstSymbolicFrontend::AtomView);
    let initial = formulation
        .initial_unknowns(TEMPERATURE_K)
        .expect("monolithic initial point");
    let expr_residual = expr
        .residual_for_test(&initial)
        .expect("ExprLegacy initial residual");
    let atom_residual = atom
        .residual_for_test(&initial)
        .expect("AtomView initial residual");
    let expr_jacobian = expr
        .jacobian_for_test(&initial)
        .expect("ExprLegacy initial Jacobian");
    let atom_jacobian = atom
        .jacobian_for_test(&initial)
        .expect("AtomView initial Jacobian");
    let (residual_delta, residual_row) =
        max_abs_vector_delta(expr_residual.as_slice(), atom_residual.as_slice());
    let (jacobian_delta, jacobian_index) = max_abs_matrix_delta(&expr_jacobian, &atom_jacobian);
    assert!(
        residual_delta <= 1.0e-12,
        "STANJAN initial residual frontend drift {residual_delta:.17e} at row {residual_row}"
    );
    assert!(
        jacobian_delta <= 1.0e-12,
        "STANJAN initial Jacobian frontend drift {jacobian_delta:.17e} at {jacobian_index:?}"
    );

    let target_enthalpy_j = request
        .constraint
        .target_enthalpy_joules()
        .expect("P,H target enthalpy")
        .joules();
    let expr_result = solve_public_ph_facade(target_enthalpy_j, RstSymbolicFrontend::ExprLegacy);
    let atom_result = solve_public_ph_facade(target_enthalpy_j, RstSymbolicFrontend::AtomView);
    assert!((expr_result.0 - TEMPERATURE_K).abs() <= 1.0e-3);
    assert!((atom_result.0 - TEMPERATURE_K).abs() <= 1.0e-3);
    assert_eq!(expr_result.2, atom_result.2, "phase topology differs");
    let expr_total = expr_result.1.iter().sum::<f64>();
    let atom_total = atom_result.1.iter().sum::<f64>();
    let fraction_delta = expr_result
        .1
        .iter()
        .zip(&atom_result.1)
        .map(|(expr, atom)| (expr / expr_total - atom / atom_total).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        fraction_delta <= 1.0e-8,
        "public facade mole-fraction frontend drift {fraction_delta:.17e}"
    );
}

#[test]
#[ignore = "diagnostic STANJAN-15 AtomView versus ExprLegacy nonlinear path; run explicitly"]
fn stanjan15_monolithic_ph_compares_symbolic_linearizations_and_iteration_traces() {
    let (request, formulation) = diagnostic_inputs();
    let expr = prepare_frontend(&formulation, RstSymbolicFrontend::ExprLegacy);
    let atom = prepare_frontend(&formulation, RstSymbolicFrontend::AtomView);
    report_tiny_constant_lambdify_probe();
    let initial = formulation
        .initial_unknowns(TEMPERATURE_K)
        .expect("monolithic seed");
    println!(
        "STANJAN-15 P,H diagnostic: x=[ln(moles), transformed_temperature], species={}, target_H={:.9e} J",
        seed_moles().len(),
        request
            .constraint
            .target_enthalpy_joules()
            .expect("P,H target")
            .joules(),
    );
    compare_at("initial_seed", &formulation, &expr, &atom, &initial);

    let expr_results = solve_with_capture(&request, RstSymbolicFrontend::ExprLegacy);
    let atom_results = solve_with_capture(&request, RstSymbolicFrontend::AtomView);
    for (frontend, result) in [
        ("ExprLegacy", expr_results.first()),
        ("AtomView", atom_results.first()),
    ] {
        if let Some((x, _)) = result {
            compare_at(frontend, &formulation, &expr, &atom, x);
        }
    }
    assert!(
        !expr_results.is_empty(),
        "ExprLegacy should reach an RST attempt"
    );
    assert!(
        !atom_results.is_empty(),
        "AtomView should reach an RST attempt"
    );
}
