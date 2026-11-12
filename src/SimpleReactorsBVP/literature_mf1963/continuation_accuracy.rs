//! Focused continuation accuracy check for the MF-HIGH hot branch.

use super::*;

#[test]
fn mf1963_damped_hot_continuation_01_to_0125_reports_independent_defect() {
    let base = cases()[1];
    let nodes = 301;
    let hot = solve_mf_correction_newton(base, nodes)
        .unwrap_or_else(|error| panic!("MF-HIGH hot control: {error}"));
    let base_seed = damped_reduced_seed_from_physical_profile(&hot, nodes);
    let base_solution = solve_mf_damped_with_seed(
        base,
        base_seed,
        nodes,
        0,
        GridRefinementMethod::DoublePoints,
        1,
        false,
        false,
    )
    .unwrap_or_else(|error| panic!("Damped delta=0.1 hot-branch baseline: {error}"));
    let baseline_defect = mf_discrete_trapezoid_defects(base, &base_solution);
    let base_temperature = base_solution
        .solver
        .get_solution()
        .expect("baseline solution")
        .column(0)
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    assert_damped_reports_converged(&base_solution);
    assert!(base_temperature >= 1000.0, "baseline left the hot branch");

    let next = Case {
        delta: 0.125,
        ..base
    };
    let next_seed = damped_reduced_seed_from_physical_profile(&base_solution, nodes);
    let continued = solve_mf_damped_with_seed(
        next,
        next_seed,
        nodes,
        0,
        GridRefinementMethod::DoublePoints,
        1,
        false,
        false,
    )
    .unwrap_or_else(|error| panic!("Damped delta=0.125 continuation: {error}"));
    let solution = continued.solver.get_solution().expect("continued solution");
    let temperature_max = solution
        .column(0)
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let defects = mf_discrete_trapezoid_defects(next, &continued);
    let profile_delta = physical_profile_difference(&continued, &base_solution);
    assert_damped_reports_converged(&continued);

    println!(
        "MF-04C Damped hot continuation: nodes={nodes} delta=0.1 -> 0.125 T_max={temperature_max:.9e} K baseline_defects=({:.6e},{:.6e},{:.6e}) continued_defects=({:.6e},{:.6e},{:.6e}) profile_delta=(dT={:.6e} K,dY={:.6e})",
        baseline_defect.0,
        baseline_defect.1,
        baseline_defect.2,
        defects.0,
        defects.1,
        defects.2,
        profile_delta.0,
        profile_delta.1,
    );

    assert!(
        temperature_max >= 1000.0,
        "continuation left the hot branch"
    );
    assert!(solution.iter().all(|value| value.is_finite()));
    assert!(defects.0.is_finite() && defects.1.is_finite() && defects.2.is_finite());
    assert!(
        defects.1 > 1e-5,
        "the previously observed p-equation defect unexpectedly passed its strict gate: {defects:?}"
    );
}
