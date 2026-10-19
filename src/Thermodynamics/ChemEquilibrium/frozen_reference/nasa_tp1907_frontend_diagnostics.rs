use super::*;

use RustedSciThe::numerical::SolverTelemetryMode;
use std::time::Instant;

fn max_abs_difference(left: &[f64], right: &[f64]) -> f64 {
    assert_eq!(left.len(), right.len());
    left.iter()
        .zip(right)
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max)
}

fn max_abs_matrix_difference(left: &nalgebra::DMatrix<f64>, right: &nalgebra::DMatrix<f64>) -> f64 {
    assert_eq!(left.shape(), right.shape());
    left.iter()
        .zip(right.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0, f64::max)
}

fn relative_to_shared_scale(difference: f64, left: &[f64], right: &[f64]) -> f64 {
    let scale = left
        .iter()
        .chain(right)
        .map(|value| value.abs())
        .fold(1.0_f64, f64::max);
    difference / scale
}

fn report_linearization(
    label: &str,
    problem: &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
    y: &[f64],
) {
    let residual = problem
        .residual_for_test(y)
        .expect("linearization residual must be finite");
    let jacobian = problem
        .jacobian_for_test(y)
        .expect("linearization Jacobian must be finite");
    let singular_values = jacobian.clone().svd(false, false).singular_values;
    let sigma_max = singular_values.iter().copied().fold(0.0_f64, f64::max);
    let sigma_min = singular_values
        .iter()
        .copied()
        .filter(|value| *value > 0.0)
        .fold(f64::INFINITY, f64::min);
    let gradient_norm = (jacobian.transpose() * residual.clone()).norm();
    println!(
        "linearization stage={label} residual_l2={:.6e} jacobian_max_abs={:.6e} sigma_max={sigma_max:.6e} sigma_min_nonzero={sigma_min:.6e} condition_estimate={:.6e} gradient_l2={gradient_norm:.6e}",
        residual.norm(),
        jacobian.iter().map(|value| value.abs()).fold(0.0, f64::max),
        sigma_max / sigma_min,
    );
}

fn compare_callbacks(
    stage: &str,
    expr: &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
    atom: &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
    y: &[f64],
) {
    let expr_residual = expr
        .residual_for_test(y)
        .unwrap_or_else(|error| panic!("ExprLegacy {stage} residual failed: {error}"));
    let atom_residual = atom
        .residual_for_test(y)
        .unwrap_or_else(|error| panic!("AtomView {stage} residual failed: {error}"));
    let expr_jacobian = expr
        .jacobian_for_test(y)
        .unwrap_or_else(|error| panic!("ExprLegacy {stage} Jacobian failed: {error}"));
    let atom_jacobian = atom
        .jacobian_for_test(y)
        .unwrap_or_else(|error| panic!("AtomView {stage} Jacobian failed: {error}"));
    let residual_abs = max_abs_difference(expr_residual.as_slice(), atom_residual.as_slice());
    let jacobian_abs = max_abs_matrix_difference(&expr_jacobian, &atom_jacobian);
    let worst_residual_row = expr_residual
        .iter()
        .zip(atom_residual.iter())
        .enumerate()
        .max_by(|(_, (left_a, right_a)), (_, (left_b, right_b))| {
            (*left_a - *right_a)
                .abs()
                .total_cmp(&(*left_b - *right_b).abs())
        })
        .map(|(index, (expr, atom))| (index, *expr, *atom));
    println!(
        "callback stage={stage} expr_residual_l2={:.6e} atom_residual_l2={:.6e} residual_max_abs={residual_abs:.6e} residual_relative={:.6e} jacobian_max_abs={jacobian_abs:.6e} worst_row={worst_residual_row:?}",
        expr_residual.norm(),
        atom_residual.norm(),
        relative_to_shared_scale(
            residual_abs,
            expr_residual.as_slice(),
            atom_residual.as_slice()
        ),
    );
}

fn solve_and_report(
    label: &str,
    solver: RustedSciTheSolver,
    problem: &crate::Thermodynamics::ChemEquilibrium::equilibrium_rst_backend::RstPreparedProblem,
    seed: &[f64],
) -> Option<(Vec<f64>, Vec<f64>)> {
    let contract = RustedSciTheSolveContract::new(1.0e-12, 500)
        .expect("diagnostic solve contract must validate")
        .with_telemetry_mode(SolverTelemetryMode::Detailed);
    let started = Instant::now();
    let result = solver.solve_with_history_for_test(problem, seed, &contract);
    match result {
        Ok(result) => {
            let stats = &result.statistics;
            println!(
                "direct_solve case={label} method={} status={:?} iterations={} residual_l2={:.6e} wall_ms={:.3} callbacks(res={},jac={}) linear(solve={},factor={}) accepted={} rejected={} residual_ms={:.3} jacobian_ms={:.3} linear_ms={:.3}",
                solver.name(),
                result.termination,
                result.iterations,
                result.residual_norm,
                started.elapsed().as_secs_f64() * 1.0e3,
                stats.residual_evaluations,
                stats.jacobian_evaluations,
                stats.linear_solves,
                stats.linear_factorizations,
                stats.accepted_steps,
                stats.rejected_steps,
                stats.residual_duration.as_secs_f64() * 1.0e3,
                stats.jacobian_duration.as_secs_f64() * 1.0e3,
                stats.linear_solve_duration.as_secs_f64() * 1.0e3,
            );
            let history_len = result.history.len();
            for (index, (history, attempt)) in
                result.history.iter().zip(stats.attempts.iter()).enumerate()
            {
                if index >= 5 && index + 5 < history_len {
                    continue;
                }
                println!(
                    "  iteration={} accepted={} residual_l2={:.6e} step_norm={:.6e} callback_calls(r={},j={}) linear(solve={},factor={}) callback_ms(r={:.3},j={:.3})",
                    history.iteration,
                    history.accepted,
                    history.residual_norm,
                    history.step_norm,
                    attempt.residual_evaluations,
                    attempt.jacobian_evaluations,
                    attempt.linear_solves,
                    attempt.linear_factorizations,
                    attempt.residual_duration.as_secs_f64() * 1.0e3,
                    attempt.jacobian_duration.as_secs_f64() * 1.0e3,
                );
            }
            Some((
                result.x.as_slice().to_vec(),
                result.residual.as_slice().to_vec(),
            ))
        }
        Err(error) => {
            println!(
                "direct_solve case={label} method={} status=ERROR error={error} wall_ms={:.3}",
                solver.name(),
                started.elapsed().as_secs_f64() * 1.0e3,
            );
            None
        }
    }
}

#[test]
fn tp1907_atom_view_failure_is_localized_by_stage_and_native_telemetry() {
    let fixture = fixture();
    let temperature_k = 700.0;
    let scale = 1.0e8;
    let physical = prepared_fixed_tp1907_problem_at_scale(&fixture, temperature_k, scale);

    let started = Instant::now();
    let (expr_physical, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
        &physical,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
        SymbolicLambdifyFrontend::ExprLegacy,
    )
    .expect("physical ExprLegacy problem must prepare");
    let expr_prepare_ms = started.elapsed().as_secs_f64() * 1.0e3;
    let started = Instant::now();
    let (atom_physical, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
        &physical,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
        SymbolicLambdifyFrontend::AtomViewNative,
    )
    .expect("physical AtomView problem must prepare");
    let atom_prepare_ms = started.elapsed().as_secs_f64() * 1.0e3;

    let initial = scaled_initial_composition(&fixture, scale);
    let (inventory, _, _) = frozen_element_inventory(&fixture, &initial);
    let options = EquilibriumSolveOptions::default()
        .with_max_iterations(500)
        .expect("diagnostic facade iteration budget must validate")
        .with_tolerance(1.0e-12)
        .expect("diagnostic facade tolerance must validate");
    let expr_outcome = with_test_rst_frontend(SymbolicLambdifyFrontend::ExprLegacy, || {
        calculator_tp1907_builder(&fixture, repository())
            .solve_options(options)
            .solver_policy(SolverPolicy::Single(SolverBackend::RustedSciThe(
                RustedSciTheSolver::TrustRegionLevenbergMarquardt,
            )))
            .expect("diagnostic facade solver policy must validate")
            .at_temperature(temperature_k)
            .element_inventory(inventory)
            .element_numerical_seed(initial.clone())
            .solve()
    })
    .expect("ExprLegacy facade route must provide the comparison solution");
    let expr_point = match expr_outcome {
        EquilibriumCalculatorOutcome::PtPoint(point) => point,
        other => panic!("expected TP-1907 diagnostic point, got {other:?}"),
    };
    let expr_y_final = expr_point
        .solution()
        .accepted_solution()
        .log_moles()
        .to_vec();

    let seed = physical.problem().initial_log_moles().as_slice();
    let midpoint = seed
        .iter()
        .zip(expr_y_final.iter())
        .map(|(&start, &end)| 0.5 * (start + end))
        .collect::<Vec<_>>();
    println!(
        "TP-1907 AtomView diagnostic scale={scale:.1e} dimensions={} preparation_ms(expr={expr_prepare_ms:.3},atom={atom_prepare_ms:.3})",
        seed.len()
    );
    compare_callbacks("physical_seed", &expr_physical, &atom_physical, seed);
    compare_callbacks(
        "physical_midpoint_to_expr_solution",
        &expr_physical,
        &atom_physical,
        &midpoint,
    );
    compare_callbacks(
        "physical_expr_facade_solution",
        &expr_physical,
        &atom_physical,
        &expr_y_final,
    );

    let _physical_expr_result = solve_and_report(
        "physical_expr",
        RustedSciTheSolver::TrustRegionLevenbergMarquardt,
        &expr_physical,
        seed,
    );
    let _physical_atom_result = solve_and_report(
        "physical_atom",
        RustedSciTheSolver::TrustRegionLevenbergMarquardt,
        &atom_physical,
        seed,
    );
    let failure_contract = RustedSciTheSolveContract::new(1.0e-12, 500)
        .expect("failure trace contract must validate")
        .with_telemetry_mode(SolverTelemetryMode::Detailed);
    match RustedSciTheSolver::TrustRegionLevenbergMarquardt.solve_with_failure_trace_for_test(
        &atom_physical,
        seed,
        &failure_contract,
    ) {
        Ok(result) => println!(
            "atom_failure_trace unexpectedly_succeeded=true termination={:?} residual_l2={:.6e}",
            result.termination, result.residual_norm
        ),
        Err(trace) => {
            let failing_y = trace
                .residual_inputs
                .last()
                .expect("failure trace must retain the failing residual input");
            let last_finite = trace
                .residual_inputs
                .iter()
                .enumerate()
                .rev()
                .find_map(|(index, candidate)| {
                    atom_physical
                        .residual_for_test(candidate)
                        .ok()
                        .map(|residual| (index, candidate.clone(), residual))
                })
                .expect("trace must contain a previous finite AtomView residual");
            compare_callbacks(
                "physical_last_finite_before_atom_failure",
                &expr_physical,
                &atom_physical,
                &last_finite.1,
            );
            report_linearization("last_finite_atom_state", &atom_physical, &last_finite.1);
            let variables = atom_physical.symbolic_payload_for_test().1;
            let species = physical.problem().species();
            let mut largest_state_coordinates = last_finite
                .1
                .iter()
                .enumerate()
                .map(|(index, value)| {
                    (
                        index,
                        variables[index].as_str(),
                        species[index].as_str(),
                        *value,
                    )
                })
                .collect::<Vec<_>>();
            largest_state_coordinates.sort_by(|left, right| right.3.abs().total_cmp(&left.3.abs()));
            let non_finite_trial_coordinates = failing_y
                .iter()
                .enumerate()
                .filter_map(|(index, value)| {
                    (!value.is_finite()).then_some((
                        index,
                        variables[index].as_str(),
                        species[index].as_str(),
                        *value,
                    ))
                })
                .collect::<Vec<_>>();
            let (failure_status, failure_values) =
                atom_physical.residual_components_for_test(failing_y);
            let non_finite_rows = failure_values
                .iter()
                .enumerate()
                .filter_map(|(index, value)| (!value.is_finite()).then_some(index))
                .collect::<Vec<_>>();
            let expr_at_failure = expr_physical.residual_for_test(failing_y);
            println!(
                "atom_failure_trace error={} residual_calls={} jacobian_calls={} failing_call_index={} last_finite_call_index={} last_finite_l2={:.6e} largest_last_finite_coordinates={:?} failing_y_all_finite={} non_finite_trial_coordinates={:?} non_finite_residual_rows={non_finite_rows:?} expr_at_same_trial={:?}",
                trace.error,
                trace.residual_inputs.len(),
                trace.jacobian_inputs.len(),
                trace.residual_inputs.len().saturating_sub(1),
                last_finite.0,
                last_finite.2.norm(),
                largest_state_coordinates.iter().take(5).collect::<Vec<_>>(),
                failing_y.iter().all(|value| value.is_finite()),
                non_finite_trial_coordinates,
                expr_at_failure.as_ref().map(|residual| residual.norm()),
            );
            println!("atom_failure_residual_status={failure_status:?}");
            let last_jacobian_point = trace
                .jacobian_inputs
                .last()
                .expect("failure trace must retain its last Jacobian evaluation");
            compare_callbacks(
                "physical_last_jacobian_before_atom_failure",
                &expr_physical,
                &atom_physical,
                last_jacobian_point,
            );
            report_linearization(
                "last_jacobian_point_expr",
                &expr_physical,
                last_jacobian_point,
            );
            report_linearization(
                "last_jacobian_point_atom",
                &atom_physical,
                last_jacobian_point,
            );
            for (index, candidate) in trace
                .residual_inputs
                .iter()
                .enumerate()
                .rev()
                .filter(|(_, candidate)| candidate.iter().all(|value| value.is_finite()))
                .take(5)
                .collect::<Vec<_>>()
                .into_iter()
                .rev()
            {
                let residual = atom_physical
                    .residual_for_test(candidate)
                    .expect("recent finite callback point must remain evaluable");
                println!(
                    "atom_residual_path call={index} l2={:.6e} max_abs_y={:.6e}",
                    residual.norm(),
                    candidate
                        .iter()
                        .map(|value| value.abs())
                        .fold(0.0, f64::max)
                );
            }
        }
    }

    let (normalized, normalization_scale) =
        normalized_fixed_problem_from_input(&fixture, temperature_k, &physical);
    let (expr_normalized, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
        &normalized,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
        SymbolicLambdifyFrontend::ExprLegacy,
    )
    .expect("normalized ExprLegacy problem must prepare");
    let (atom_normalized, _) = prepare_rst_symbolic_problem_from_prepared_with_frontend(
        &normalized,
        crate::Thermodynamics::ChemEquilibrium::equilibrium_timing::EquilibriumTimingMode::Disabled,
        SymbolicLambdifyFrontend::AtomViewNative,
    )
    .expect("normalized AtomView problem must prepare");
    let normalized_seed = normalized.problem().initial_log_moles().as_slice();
    println!("normalization physical_inventory_scale={normalization_scale:.9e}");
    compare_callbacks(
        "normalized_seed",
        &expr_normalized,
        &atom_normalized,
        normalized_seed,
    );
    let normalized_expr_result = solve_and_report(
        "normalized_expr",
        RustedSciTheSolver::TrustRegionLevenbergMarquardt,
        &expr_normalized,
        normalized_seed,
    );
    let normalized_atom_result = solve_and_report(
        "normalized_atom",
        RustedSciTheSolver::TrustRegionLevenbergMarquardt,
        &atom_normalized,
        normalized_seed,
    );
    if let (Some((expr_y, expr_residual)), Some((atom_y, atom_residual))) =
        (normalized_expr_result, normalized_atom_result)
    {
        let (_, variables, _, _) = expr_normalized.symbolic_payload_for_test();
        let species = normalized.problem().species();
        let mut coordinate_deltas = expr_y
            .iter()
            .zip(atom_y.iter())
            .zip(variables.iter())
            .zip(species.iter())
            .enumerate()
            .map(|(index, (((expr, atom), variable), species))| {
                (
                    index,
                    variable.as_str(),
                    species.as_str(),
                    (expr - atom).abs(),
                    *expr,
                    *atom,
                )
            })
            .collect::<Vec<_>>();
        coordinate_deltas.sort_by(|left, right| right.3.total_cmp(&left.3));
        println!(
            "normalized_terminal_comparison max_abs_y_delta={:.6e} max_abs_residual_delta={:.6e} largest_coordinate_deltas={:?}",
            max_abs_difference(&expr_y, &atom_y),
            max_abs_difference(&expr_residual, &atom_residual),
            coordinate_deltas.iter().take(5).collect::<Vec<_>>(),
        );
        compare_callbacks(
            "normalized_expr_terminal_y",
            &expr_normalized,
            &atom_normalized,
            &expr_y,
        );
        compare_callbacks(
            "normalized_atom_terminal_y",
            &expr_normalized,
            &atom_normalized,
            &atom_y,
        );
    }
}
