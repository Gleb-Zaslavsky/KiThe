#[cfg(test)]
mod tests {
    use nalgebra::{DMatrix, DVector};
    use std::fs;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_activity::PhaseActivityModel;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::{Phase, R};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_phase_stability::{
        IdealTpdProblem, fit_element_potentials,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        InitialPhaseSet, PhaseManager, PhaseSet, PhaseStatus, PhaseTransitionPlan,
        compute_phase_stability_reports, compute_phase_totals,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_ternary_vle_antoine_gauge::{
        FROZEN_ANTOINE_REFERENCE_PRESSURE_PA, FrozenAntoineGaugeState,
        TernaryRaoultFlashPhaseClass, common_temperature_interval, gauge_standard_gibbs_j_mol,
        load_nist_ternary_vle_antoine_gauge, load_nist_ternary_vle_interior_rows,
        molecular_element_matrix, phase_transfer_reaction_basis, saturation_pressure_from_gauge_pa,
        saturation_pressure_pa, solve_raoult_bubble_temperature, solve_raoult_flash,
        ternary_gauge_gibbs_functions,
    };
    use crate::Thermodynamics::ChemEquilibrium::prelude::PhaseIndex;

    #[test]
    fn frozen_antoine_universe_has_three_independent_molecules_and_three_transfer_directions() {
        let dataset = load_nist_ternary_vle_antoine_gauge().unwrap();
        assert!(dataset.is_external_evidence());
        let (labels, molecular) = molecular_element_matrix(dataset.rows()).unwrap();
        assert_eq!(labels, ["C", "Cl", "H"]);
        assert_eq!(molecular.nrows(), 3);
        assert_eq!(molecular.ncols(), 3);
        let (rank, phase_transfer_reactions) =
            phase_transfer_reaction_basis(dataset.rows()).unwrap();
        assert_eq!(rank, 3);
        assert_eq!(phase_transfer_reactions, 3);
    }

    #[test]
    fn common_interval_is_derived_from_all_three_frozen_antoine_records() {
        let dataset = load_nist_ternary_vle_antoine_gauge().unwrap();
        assert_eq!(
            common_temperature_interval(dataset.rows()).unwrap(),
            (335.19, 384.66)
        );
    }

    #[test]
    fn pure_component_gibbs_gauge_round_trips_each_antoine_pressure_without_vle_data() {
        let dataset = load_nist_ternary_vle_antoine_gauge().unwrap();
        for record in dataset.rows() {
            for temperature in [340.0, 350.0, 360.0, 370.0, 380.0] {
                let pressure = saturation_pressure_pa(record, temperature).unwrap();
                let gas =
                    gauge_standard_gibbs_j_mol(record, FrozenAntoineGaugeState::Gas, temperature)
                        .unwrap();
                let liquid = gauge_standard_gibbs_j_mol(
                    record,
                    FrozenAntoineGaugeState::Liquid,
                    temperature,
                )
                .unwrap();
                let recovered =
                    saturation_pressure_from_gauge_pa(gas, liquid, temperature).unwrap();
                assert!(
                    (recovered - pressure).abs() <= pressure * 1.0e-13,
                    "{} at {temperature} K",
                    record.compound
                );
                assert!(pressure < FROZEN_ANTOINE_REFERENCE_PRESSURE_PA);
                assert!(
                    liquid < gas,
                    "{} below its normal boiling point",
                    record.compound
                );
            }
        }
    }

    #[test]
    fn frozen_antoine_gauge_rejects_extrapolation_and_never_rewrites_evidence() {
        let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("nist_webbook");
        let metadata = directory.join("toluene_ethylbenzene_chlorobenzene_antoine.metadata.json");
        let rows = directory.join("toluene_ethylbenzene_chlorobenzene_antoine.rows.json");
        let thermoml_directory = directory.join("..").join("nist_thermoml");
        let thermoml_metadata =
            thermoml_directory.join("toluene_ethylbenzene_chlorobenzene_selected.metadata.json");
        let thermoml_rows =
            thermoml_directory.join("toluene_ethylbenzene_chlorobenzene_selected.rows.json");
        let before_metadata = fs::read(&metadata).unwrap();
        let before_rows = fs::read(&rows).unwrap();
        let before_thermoml_metadata = fs::read(&thermoml_metadata).unwrap();
        let before_thermoml_rows = fs::read(&thermoml_rows).unwrap();
        let dataset = load_nist_ternary_vle_antoine_gauge().unwrap();
        let _ = load_nist_ternary_vle_interior_rows().unwrap();

        assert!(saturation_pressure_pa(&dataset.rows()[0], 300.0).is_err());
        assert!(saturation_pressure_pa(&dataset.rows()[2], 410.0).is_err());
        assert_eq!(fs::read(metadata).unwrap(), before_metadata);
        assert_eq!(fs::read(rows).unwrap(), before_rows);
        assert_eq!(
            fs::read(thermoml_metadata).unwrap(),
            before_thermoml_metadata
        );
        assert_eq!(fs::read(thermoml_rows).unwrap(), before_thermoml_rows);
    }

    #[test]
    fn selected_thermoml_rows_are_interior_evidence_and_raoult_derives_a_vapor_simplex() {
        let records = load_nist_ternary_vle_antoine_gauge().unwrap();
        let rows = load_nist_ternary_vle_interior_rows().unwrap();
        let interval = common_temperature_interval(records.rows()).unwrap();

        assert!(rows.is_external_evidence());
        assert_eq!(rows.rows().len(), 4);
        for row in rows.rows() {
            assert!(
                row.source_temperature_k >= interval.0 && row.source_temperature_k <= interval.1,
                "source row at {} K must stay in the common Antoine interval",
                row.source_temperature_k
            );
            let solution = solve_raoult_bubble_temperature(
                records.rows(),
                [
                    row.liquid_toluene_mole_fraction,
                    row.liquid_ethylbenzene_mole_fraction,
                    row.liquid_chlorobenzene_mole_fraction,
                ],
                row.pressure_pa,
            )
            .unwrap();
            assert!(solution.temperature_k >= interval.0 && solution.temperature_k <= interval.1);
            assert_eq!(solution.pressure_pa, row.pressure_pa);
            assert!(solution.iterations > 0);
            assert!(
                solution
                    .vapor_mole_fractions
                    .iter()
                    .all(|value| value.is_finite() && *value > 0.0 && *value < 1.0)
            );
            assert!((solution.vapor_mole_fractions.iter().sum::<f64>() - 1.0).abs() <= 1.0e-10);
        }
    }

    #[test]
    fn canonical_tpd_recovers_both_three_component_boundary_minima_from_the_antoine_gauge() {
        let records = load_nist_ternary_vle_antoine_gauge().unwrap();
        let rows = load_nist_ternary_vle_interior_rows().unwrap();
        let (_, element_matrix) = molecular_element_matrix(records.rows()).unwrap();

        for row in rows.rows() {
            let liquid = [
                row.liquid_toluene_mole_fraction,
                row.liquid_ethylbenzene_mole_fraction,
                row.liquid_chlorobenzene_mole_fraction,
            ];
            let boundary =
                solve_raoult_bubble_temperature(records.rows(), liquid, row.pressure_pa).unwrap();
            let temperature = boundary.temperature_k;
            let rt = R * temperature;
            let gas_chemical_potentials = boundary
                .vapor_mole_fractions
                .iter()
                .map(|&mole_fraction| {
                    rt * (mole_fraction * row.pressure_pa / FROZEN_ANTOINE_REFERENCE_PRESSURE_PA)
                        .ln()
                })
                .collect::<Vec<_>>();
            let gas_fit = fit_element_potentials(
                &element_matrix,
                &DVector::from_vec(gas_chemical_potentials),
                &[0, 1, 2],
            )
            .unwrap();
            assert!(gas_fit.max_abs_residual() <= gas_fit.residual_tolerance());

            let liquid_standard_gibbs = records
                .rows()
                .iter()
                .map(|record| {
                    gauge_standard_gibbs_j_mol(record, FrozenAntoineGaugeState::Liquid, temperature)
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let liquid_candidate = IdealTpdProblem::new(
                PhaseActivityModel::IdealSolution,
                liquid_standard_gibbs.clone(),
                element_matrix.clone(),
                gas_fit.potentials().to_vec(),
                temperature,
                row.pressure_pa,
                FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
            )
            .unwrap();
            let liquid_minimum = liquid_candidate
                .constrained_minimum(&element_matrix)
                .unwrap();
            assert!(liquid_minimum.minimum().minimum_tpd().abs() <= 1.0e-7);
            assert_eq!(liquid_minimum.active_component_count(), 3);
            assert_eq!(liquid_minimum.independent_constraint_count(), 0);
            for (actual, expected) in liquid_minimum
                .minimum()
                .incipient_composition()
                .iter()
                .zip(liquid)
            {
                assert!((actual - expected).abs() <= 1.0e-8);
            }

            let liquid_chemical_potentials = liquid_standard_gibbs
                .iter()
                .zip(liquid)
                .map(|(&standard_gibbs, mole_fraction)| standard_gibbs + rt * mole_fraction.ln())
                .collect::<Vec<_>>();
            let liquid_fit = fit_element_potentials(
                &element_matrix,
                &DVector::from_vec(liquid_chemical_potentials),
                &[0, 1, 2],
            )
            .unwrap();
            assert!(liquid_fit.max_abs_residual() <= liquid_fit.residual_tolerance());
            let gas_candidate = IdealTpdProblem::new(
                PhaseActivityModel::IdealGas,
                vec![0.0; 3],
                element_matrix.clone(),
                liquid_fit.potentials().to_vec(),
                temperature,
                row.pressure_pa,
                FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
            )
            .unwrap();
            let gas_minimum = gas_candidate.constrained_minimum(&element_matrix).unwrap();
            assert!(gas_minimum.minimum().minimum_tpd().abs() <= 1.0e-7);
            assert_eq!(gas_minimum.active_component_count(), 3);
            assert_eq!(gas_minimum.independent_constraint_count(), 0);
            for (actual, expected) in gas_minimum
                .minimum()
                .incipient_composition()
                .iter()
                .zip(boundary.vapor_mole_fractions)
            {
                assert!((actual - expected).abs() <= 1.0e-8);
            }
        }
    }

    #[test]
    fn canonical_phase_control_decision_tracks_both_ternary_boundary_directions() {
        let records = load_nist_ternary_vle_antoine_gauge().unwrap();
        let row = load_nist_ternary_vle_interior_rows().unwrap().rows()[2];
        let liquid = [
            row.liquid_toluene_mole_fraction,
            row.liquid_ethylbenzene_mole_fraction,
            row.liquid_chlorobenzene_mole_fraction,
        ];
        let boundary =
            solve_raoult_bubble_temperature(records.rows(), liquid, row.pressure_pa).unwrap();
        let (_, molecular_elements) = molecular_element_matrix(records.rows()).unwrap();
        let elements = DMatrix::from_fn(6, molecular_elements.ncols(), |row, column| {
            molecular_elements[(row % 3, column)]
        });
        let gibbs = ternary_gauge_gibbs_functions(records.rows()).unwrap();
        let phases = vec![
            Phase {
                kind: PhaseActivityModel::IdealGas,
                species: vec![0, 1, 2],
            },
            Phase {
                kind: PhaseActivityModel::IdealSolution,
                species: vec![3, 4, 5],
            },
        ];
        let species_phase = [0, 0, 0, 1, 1, 1];

        for (gas_is_active, pressure_factor, expected_activation, expected_composition) in [
            (true, 0.98, None, None),
            (true, 1.02, Some(1usize), Some(liquid)),
            (
                false,
                0.98,
                Some(0usize),
                Some(boundary.vapor_mole_fractions),
            ),
            (false, 1.02, None, None),
        ] {
            let pressure = row.pressure_pa * pressure_factor;
            let active_composition = if gas_is_active {
                boundary.vapor_mole_fractions
            } else {
                liquid
            };
            let log_moles = (0..6)
                .map(|component| {
                    let is_active_component = (component < 3) == gas_is_active;
                    if is_active_component {
                        active_composition[component % 3].ln()
                    } else {
                        1.0e-30_f64.ln()
                    }
                })
                .collect::<Vec<_>>();
            let active_phase = usize::from(!gas_is_active);
            let phase_set = PhaseSet::from_policy(
                &InitialPhaseSet::Explicit {
                    active: vec![PhaseIndex::new(active_phase, 2).unwrap()],
                    excluded: Vec::new(),
                },
                &[gas_is_active, !gas_is_active],
            )
            .unwrap();
            let reports = compute_phase_stability_reports(
                &log_moles,
                &gibbs,
                &phases,
                &species_phase,
                &elements,
                boundary.temperature_k,
                pressure,
                FROZEN_ANTOINE_REFERENCE_PRESSURE_PA,
                &phase_set,
            )
            .unwrap();
            let candidate_phase = 1 - active_phase;
            let candidate = &reports[candidate_phase];
            let minimum_tpd = candidate.minimum_tpd.unwrap();
            let plan = PhaseManager::default()
                .classify_phases_at_temperature(
                    boundary.temperature_k,
                    &compute_phase_totals(&log_moles, &species_phase),
                    &reports,
                    &phase_set,
                )
                .unwrap();

            match expected_activation {
                Some(phase) => {
                    assert!(minimum_tpd < 0.0, "pressure factor {pressure_factor}");
                    assert_eq!(
                        plan,
                        PhaseTransitionPlan::Activate {
                            phase: PhaseIndex::new(phase, 2).unwrap()
                        }
                    );
                    let incipient = candidate.incipient_composition.as_ref().unwrap();
                    assert_eq!(incipient.len(), 3);
                    for (actual, expected) in incipient.iter().zip(expected_composition.unwrap()) {
                        assert!((actual - expected).abs() <= 1.0e-8);
                    }
                    let mut transitioned = phase_set.clone();
                    transitioned.activate(PhaseIndex::new(phase, 2).unwrap());
                    assert_eq!(
                        transitioned.status(PhaseIndex::new(phase, 2).unwrap()),
                        PhaseStatus::Appeared
                    );
                }
                None => {
                    assert!(minimum_tpd > 0.0, "pressure factor {pressure_factor}");
                    assert!(matches!(plan, PhaseTransitionPlan::NoTransition { .. }));
                }
            }
        }
    }

    #[test]
    fn independent_ternary_rachford_rice_flash_reconstructs_the_interior_split_and_endpoints() {
        let records = load_nist_ternary_vle_antoine_gauge().unwrap();
        let row = load_nist_ternary_vle_interior_rows().unwrap().rows()[3];
        let liquid = [
            row.liquid_toluene_mole_fraction,
            row.liquid_ethylbenzene_mole_fraction,
            row.liquid_chlorobenzene_mole_fraction,
        ];
        let boundary =
            solve_raoult_bubble_temperature(records.rows(), liquid, row.pressure_pa).unwrap();
        let expected_vapor_fraction = 0.5;
        let bulk = std::array::from_fn(|index| {
            expected_vapor_fraction * boundary.vapor_mole_fractions[index]
                + (1.0 - expected_vapor_fraction) * liquid[index]
        });

        let two_phase = solve_raoult_flash(
            records.rows(),
            boundary.temperature_k,
            row.pressure_pa,
            bulk,
        )
        .unwrap();
        assert_eq!(
            two_phase.phase_class,
            TernaryRaoultFlashPhaseClass::TwoPhase
        );
        assert!(two_phase.iterations > 0);
        assert!((two_phase.vapor_fraction.unwrap() - expected_vapor_fraction).abs() <= 1.0e-10);
        for (actual, expected) in two_phase.liquid_mole_fractions.unwrap().iter().zip(liquid) {
            assert!((actual - expected).abs() <= 1.0e-10);
        }
        for (actual, expected) in two_phase
            .vapor_mole_fractions
            .unwrap()
            .iter()
            .zip(boundary.vapor_mole_fractions)
        {
            assert!((actual - expected).abs() <= 1.0e-10);
        }

        let all_liquid = solve_raoult_flash(
            records.rows(),
            boundary.temperature_k,
            row.pressure_pa * 2.0,
            bulk,
        )
        .unwrap();
        assert_eq!(
            all_liquid.phase_class,
            TernaryRaoultFlashPhaseClass::AllLiquid
        );
        assert_eq!(all_liquid.liquid_mole_fractions, Some(bulk));
        assert_eq!(all_liquid.vapor_mole_fractions, None);

        let all_vapor = solve_raoult_flash(
            records.rows(),
            boundary.temperature_k,
            row.pressure_pa * 0.5,
            bulk,
        )
        .unwrap();
        assert_eq!(
            all_vapor.phase_class,
            TernaryRaoultFlashPhaseClass::AllVapor
        );
        assert_eq!(all_vapor.liquid_mole_fractions, None);
        assert_eq!(all_vapor.vapor_mole_fractions, Some(bulk));
    }

    /// External characterization only: NIST ThermoML was not used to tune the
    /// frozen Antoine records, so this comparison reports rather than accepts
    /// a source/model difference.
    #[test]
    #[ignore = "release-oriented NIST ternary VLE versus independent Antoine/Raoult characterization"]
    fn i5_nist_ternary_vle_antoine_raoult_temperature_characterization() {
        let records = load_nist_ternary_vle_antoine_gauge().unwrap();
        let rows = load_nist_ternary_vle_interior_rows().unwrap();
        let interval = common_temperature_interval(records.rows()).unwrap();
        let mut squared_error_sum = 0.0;
        let mut max_absolute_error = 0.0_f64;

        println!("NIST ThermoML ternary VLE Antoine/Raoult temperature characterization");
        println!("source=T/P/x only; vapor composition is derived, never read from ThermoML");
        println!(
            "  x_tol   x_ethyl   x_chl      P kPa   T_NIST K  T_Raoult K  delta K  relative %"
        );
        for row in rows.rows() {
            let solution = solve_raoult_bubble_temperature(
                records.rows(),
                [
                    row.liquid_toluene_mole_fraction,
                    row.liquid_ethylbenzene_mole_fraction,
                    row.liquid_chlorobenzene_mole_fraction,
                ],
                row.pressure_pa,
            )
            .unwrap();
            let difference = solution.temperature_k - row.source_temperature_k;
            let relative_percent = 100.0 * difference / row.source_temperature_k;
            assert!(solution.temperature_k >= interval.0 && solution.temperature_k <= interval.1);
            assert!(difference.is_finite() && relative_percent.is_finite());
            squared_error_sum += difference * difference;
            max_absolute_error = max_absolute_error.max(difference.abs());
            println!(
                "  {:>5.3}   {:>5.3}   {:>5.3}   {:>8.3}   {:>8.3}   {:>10.3}  {:+7.3}  {:+9.4}",
                row.liquid_toluene_mole_fraction,
                row.liquid_ethylbenzene_mole_fraction,
                row.liquid_chlorobenzene_mole_fraction,
                row.pressure_pa / 1_000.0,
                row.source_temperature_k,
                solution.temperature_k,
                difference,
                relative_percent,
            );
        }
        let rms_error = (squared_error_sum / rows.rows().len() as f64).sqrt();
        assert!(
            rms_error < 1.5,
            "frozen ThermoML characterization drifted beyond its conservative RMS guard: {rms_error} K"
        );
        assert!(
            max_absolute_error < 2.0,
            "frozen ThermoML characterization drifted beyond its conservative maximum-error guard: {max_absolute_error} K"
        );
        println!(
            "summary: rows={} RMS delta={rms_error:.4} K max abs delta={max_absolute_error:.4} K",
            rows.rows().len()
        );
    }
}
