//! Feasibility diagnostics for transactional rollback in real continuation.
//!
//! The first matrix is a natural-failure feasibility diagnostic. The second
//! is a real-data rollback story using a narrowly scoped test-only transition
//! boundary failpoint after that feasibility search found no natural failure.

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_range::{
        PhEnthalpyGrid, PhRangeError, PhRangeRequest, PhRangeSolution,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::PhSolveMode;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1906_chon_graphite_ph::ResolvedNasaTp1906ChonGraphitePhFixture;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::EquilibriumSolveOptions;
    use crate::Thermodynamics::ChemEquilibrium::prelude::PhaseControlPolicy;
    use crate::Thermodynamics::ChemEquilibrium::prepared_phase_control_runner::PreparedPhaseControlRunner;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};

    fn repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository().expect("bundled local repository must be available")
    }

    fn range_request<'a>(
        fixture: &'a ResolvedNasaTp1906ChonGraphitePhFixture,
        targets: &PhEnthalpyGrid,
    ) -> PhRangeRequest<'a> {
        PhRangeRequest::from_resolved_thermochemistry(
            fixture.tp1907().resolved(),
            fixture
                .tp1907()
                .initial_composition()
                .expect("initial composition must resolve"),
            101_325.0,
            101_325.0,
            targets.clone(),
            fixture.tp1907().thermochemistry().temperature_bounds(),
            900.0,
            fixture.tp1907().thermochemistry().clone(),
        )
        .expect("real P,H range request must validate")
        .with_solve_options(
            EquilibriumSolveOptions::default()
                .with_max_iterations(100)
                .expect("rollback story iteration budget must validate"),
        )
        .with_phase_control_policy(PhaseControlPolicy::default())
        .with_ph_solve_mode(PhSolveMode::NestedTemperature)
    }

    fn assert_same_physical_range(left: &PhRangeSolution, right: &PhRangeSolution) {
        assert_eq!(left.points().len(), right.points().len());
        for (left_point, right_point) in left.points().iter().zip(right.points()) {
            assert_eq!(
                left_point.report().target_enthalpy_joules().to_bits(),
                right_point.report().target_enthalpy_joules().to_bits()
            );
            assert!(
                (left_point.solution().temperature() - right_point.solution().temperature()).abs()
                    <= 1.0e-6,
                "rollback retry changed accepted temperature"
            );
            for (&left_moles, &right_moles) in left_point
                .solution()
                .equilibrium()
                .component_moles()
                .iter()
                .zip(right_point.solution().equilibrium().component_moles())
            {
                let scale = left_moles.abs().max(right_moles.abs()).max(1.0e-12);
                assert!(
                    (left_moles - right_moles).abs() / scale <= 1.0e-6,
                    "rollback retry changed accepted physical composition"
                );
            }
        }
    }

    /// Search only real TP-1906/1907 continuation points. The output is
    /// evidence for choosing a later rollback regression, not an acceptance
    /// contract: no production failure is manufactured here.
    #[test]
    #[ignore = "release feasibility matrix for real multicomponent continuation rollback"]
    fn i5_tp1906_multicomponent_continuation_rollback_feasibility_matrix() {
        let fixture = ResolvedNasaTp1906ChonGraphitePhFixture::resolve_offline(repository())
            .expect("reviewed TP-1906/1907 fixture must resolve offline");
        let cases = fixture.cases().expect("TP-1906/1907 rows must join");
        let selected = [700.0, 720.0, 740.0]
            .into_iter()
            .map(|temperature| {
                cases
                    .iter()
                    .find(|case| (case.enthalpy.temperature_k - temperature).abs() < 1e-12)
                    .expect("diagnostic must retain the requested source row")
            })
            .collect::<Vec<_>>();
        let targets =
            PhEnthalpyGrid::new(selected.iter().map(|case| case.target_enthalpy_j).collect())
                .expect("diagnostic target grid must be monotone");

        println!("NASA TP-1906/1907 continuation rollback feasibility");
        println!("budget  status       failing_point  accepted_points  detail");
        let mut found_post_acceptance_failure = false;
        let mut outcomes = Vec::new();
        for budget in [1_usize, 2, 4, 8, 16, 32, 64, 100] {
            let request = PhRangeRequest::from_resolved_thermochemistry(
                fixture.tp1907().resolved(),
                fixture
                    .tp1907()
                    .initial_composition()
                    .expect("initial composition must resolve"),
                101_325.0,
                101_325.0,
                targets.clone(),
                fixture.tp1907().thermochemistry().temperature_bounds(),
                900.0,
                fixture.tp1907().thermochemistry().clone(),
            )
            .expect("real P,H range request must validate")
            .with_solve_options(
                EquilibriumSolveOptions::default()
                    .with_max_iterations(budget)
                    .expect("diagnostic iteration budget must validate"),
            )
            .with_phase_control_policy(PhaseControlPolicy::default())
            .with_ph_solve_mode(PhSolveMode::NestedTemperature);

            match request.solve() {
                Ok(range) => {
                    assert_eq!(range.points().len(), 3);
                    outcomes.push((budget, None));
                    println!(
                        "{budget:>6}  OK           -               {:>14}  complete",
                        range.points().len()
                    );
                }
                Err(PhRangeError::Point(error)) => {
                    let index = error.index();
                    let accepted = index;
                    let post_acceptance = index > 0;
                    found_post_acceptance_failure |= post_acceptance;
                    outcomes.push((budget, Some(index)));
                    println!(
                        "{budget:>6}  FAILED       {index:>13}  {accepted:>14}  post_acceptance={post_acceptance}"
                    );
                }
                Err(error) => {
                    panic!("feasibility matrix must fail at a typed range point: {error}")
                }
            }
        }

        println!("natural post-acceptance failure found: {found_post_acceptance_failure}");
        assert_eq!(
            outcomes,
            vec![
                (1, Some(0)),
                (2, Some(0)),
                (4, Some(0)),
                (8, Some(0)),
                (16, None),
                (32, None),
                (64, None),
                (100, None),
            ],
            "the recorded feasibility boundary is a regression contract, not merely output"
        );
        assert!(
            !found_post_acceptance_failure,
            "a newly discovered natural post-acceptance failure needs a separate rollback diagnosis"
        );
    }

    /// Real-data A/B/C evidence for a failed continuation point.
    ///
    /// Route A is a clean range. Route B injects one failure after all phase
    /// transitions of the accepted first point, so the second point fails; a
    /// fresh retry then acts as the recovered route. Route C is an independent
    /// clean range. The comparison is deliberately physical and ignores
    /// diagnostic counters and the rejected attempt's internal history.
    #[test]
    #[ignore = "release real-data continuation rollback story"]
    fn i5_tp1906_multicomponent_continuation_rollback_real_route_matrix() {
        let fixture = ResolvedNasaTp1906ChonGraphitePhFixture::resolve_offline(repository())
            .expect("reviewed TP-1906/1907 fixture must resolve offline");
        let cases = fixture.cases().expect("TP-1906/1907 rows must join");
        let selected = [700.0, 720.0, 740.0]
            .into_iter()
            .map(|temperature| {
                cases
                    .iter()
                    .find(|case| (case.enthalpy.temperature_k - temperature).abs() < 1.0e-12)
                    .expect("rollback story must retain each requested source row")
            })
            .collect::<Vec<_>>();
        let targets =
            PhEnthalpyGrid::new(selected.iter().map(|case| case.target_enthalpy_j).collect())
                .expect("rollback story targets must be monotone");

        let clean_a = range_request(&fixture, &targets)
            .solve()
            .expect("clean route A must solve");
        let first_transitions = clean_a.points()[0].report().phase_control_transitions();
        assert!(
            first_transitions > 0,
            "the first real continuation point must exercise phase control"
        );

        PreparedPhaseControlRunner::arm_test_failpoint_after_transitions(first_transitions);
        let failed = range_request(&fixture, &targets).solve();
        let error = match failed {
            Err(PhRangeError::Point(error)) => error,
            other => panic!("expected a post-acceptance point failure, got {other:?}"),
        };
        assert_eq!(error.index(), 1, "failure must occur after point S0");
        assert_eq!(error.source().kind(), crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentErrorKind::AllBackendsFailed);

        let recovered_b = range_request(&fixture, &targets)
            .solve()
            .expect("one-shot failure must not poison the retry");
        let clean_c = range_request(&fixture, &targets)
            .solve()
            .expect("fresh route C must solve");
        assert_same_physical_range(&clean_a, &recovered_b);
        assert_same_physical_range(&clean_a, &clean_c);

        println!("NASA TP-1906/1907 real continuation rollback route matrix");
        println!("route  status                     points  first_transitions  failure_point");
        println!(
            "A      CLEAN                     {:>6}  {:>17}  -",
            clean_a.points().len(),
            first_transitions
        );
        println!(
            "B      FAILED@S1 -> RECOVERED   {:>6}  {:>17}  {}",
            recovered_b.points().len(),
            first_transitions,
            error.index()
        );
        println!(
            "C      FRESH CLEAN              {:>6}  {:>17}  -",
            clean_c.points().len(),
            first_transitions
        );
        println!("verdict: TransactionalRollbackInvariant");
    }
}
