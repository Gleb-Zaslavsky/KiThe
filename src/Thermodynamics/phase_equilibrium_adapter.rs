//! Pure adapters from resolved phase thermodynamics to equilibrium equations.
//!
//! The functions here do not perform database lookup and do not mutate phase
//! caches. They consume a `SystemLayout`, a borrowed property view, and an
//! element matrix whose rows are aligned with the layout components.

use crate::Thermodynamics::User_PhaseOrSolution::{PhaseLagrangeFunction, PhaseThermoPropertyView};
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DMatrix;
use std::collections::HashMap;

const R: f64 = 8.314;

/// Legacy callable shape retained at the solver boundary while the phase API
/// transitions to borrowed typed compositions.
pub type LegacyGibbsFunction = Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>;
pub type LegacyPhaseGibbsFunctions = HashMap<Option<String>, HashMap<String, LegacyGibbsFunction>>;

/// Builds the symbolic stationarity equations for Gibbs-energy minimisation.
///
/// `element_matrix` has shape `(component_count, element_count)` and its rows
/// must follow `layout.components()` exactly. This alignment is what lets a
/// substance such as water occur independently in gas and liquid phases.
pub fn build_symbolic_lagrange_equations(
    layout: &SystemLayout,
    gibbs: PhaseThermoPropertyView<'_, Expr>,
    element_matrix: &DMatrix<f64>,
    reference_temperature: f64,
) -> SubsDataResult<Vec<Expr>> {
    if !reference_temperature.is_finite() || reference_temperature <= 0.0 {
        return Err(SubsDataError::calculation_failed(
            "multiple",
            "Lagrange equations",
            format!("invalid reference temperature {}", reference_temperature),
        ));
    }
    if element_matrix.nrows() != layout.component_count() {
        return Err(SubsDataError::calculation_failed(
            "multiple",
            "Lagrange equations",
            format!(
                "element matrix has {} component rows, but the layout has {} components",
                element_matrix.nrows(),
                layout.component_count()
            ),
        ));
    }

    let multipliers = Expr::IndexedVars(element_matrix.ncols(), "Lambda").0;
    let temperature = Expr::Const(reference_temperature);
    let gas_constant = Expr::Const(R);
    let mut equations = Vec::with_capacity(layout.component_count());

    for (component_index, component) in layout.components().iter().enumerate() {
        let phase_key = component.phase.as_option();
        let phase_gibbs = gibbs
            .get(phase_key)
            .ok_or_else(|| SubsDataError::MissingData {
                field: "phase Gibbs data".to_string(),
                substance: phase_key.clone().unwrap_or_else(|| "None".to_string()),
            })?;
        let component_gibbs =
            phase_gibbs
                .get(&component.substance)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "Gibbs expression".to_string(),
                    substance: component.label(),
                })?;
        let element_term = (0..element_matrix.ncols())
            .map(|element_index| {
                Expr::Const(element_matrix[(component_index, element_index)])
                    * multipliers[element_index].clone()
            })
            .fold(Expr::Const(0.0), |sum, term| sum + term)
            .simplify();
        equations.push(
            (element_term + component_gibbs.clone() / (gas_constant.clone() * temperature.clone()))
                .simplify(),
        );
    }
    Ok(equations)
}

/// Builds numeric stationarity equations from the canonical global component
/// layout. Every component uses its own global matrix row, even when two
/// phases contain the same substance or have unequal component counts.
///
/// The public closure retains the legacy solver signature for now. Its input
/// vector is nevertheless validated against `layout.component_count()` before
/// any Gibbs callable is invoked.
pub fn build_numeric_lagrange_equations(
    layout: &SystemLayout,
    gibbs: LegacyPhaseGibbsFunctions,
    element_matrix: DMatrix<f64>,
    reference_temperature: f64,
) -> SubsDataResult<PhaseLagrangeFunction> {
    if !reference_temperature.is_finite() || reference_temperature <= 0.0 {
        return Err(SubsDataError::calculation_failed(
            "phase system",
            "numeric Lagrange equations",
            format!("invalid reference temperature {}", reference_temperature),
        ));
    }
    if element_matrix.nrows() != layout.component_count() {
        return Err(SubsDataError::calculation_failed(
            "phase system",
            "numeric Lagrange equations",
            format!(
                "element matrix has {} component rows, but the layout has {} components",
                element_matrix.nrows(),
                layout.component_count()
            ),
        ));
    }
    for component in layout.components() {
        let phase = component.phase.as_option();
        let phase_gibbs = gibbs.get(phase).ok_or_else(|| SubsDataError::MissingData {
            field: "phase Gibbs functions".to_string(),
            substance: phase.clone().unwrap_or_else(|| "None".to_string()),
        })?;
        if !phase_gibbs.contains_key(&component.substance) {
            return Err(SubsDataError::MissingData {
                field: "Gibbs function".to_string(),
                substance: component.label(),
            });
        }
    }

    let multiplier_count = element_matrix.ncols();
    let layout = layout.clone();
    Ok(Box::new(
        move |temperature, component_amounts, phase_total, multipliers| {
            if !temperature.is_finite()
                || temperature <= 0.0
                || multipliers.len() != multiplier_count
                || component_amounts
                    .as_ref()
                    .is_some_and(|amounts| amounts.len() != layout.component_count())
            {
                return Vec::new();
            }

            layout
                .components()
                .iter()
                .enumerate()
                .map(|(component_index, component)| {
                    let phase = component.phase.as_option();
                    let Some(gibbs_function) = gibbs
                        .get(phase)
                        .and_then(|phase_gibbs| phase_gibbs.get(&component.substance))
                    else {
                        return f64::NAN;
                    };
                    let constraint_term = (0..multiplier_count)
                        .map(|element_index| {
                            element_matrix[(component_index, element_index)]
                                * multipliers[element_index]
                        })
                        .sum::<f64>();
                    constraint_term
                        + gibbs_function(temperature, component_amounts.clone(), phase_total)
                            / (R * reference_temperature)
                })
                .collect()
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::User_PhaseOrSolution::PhaseSystem;
    use crate::Thermodynamics::User_substances::SubsData;
    use std::collections::HashMap;

    #[test]
    fn symbolic_adapter_uses_global_component_indices_across_phases() {
        let mut system = PhaseSystem::new();
        system.replace_dG_sym(HashMap::from([
            (
                Some("gas".to_string()),
                HashMap::from([("H2O".to_string(), Expr::Const(1.0))]),
            ),
            (
                Some("liquid".to_string()),
                HashMap::from([("H2O".to_string(), Expr::Const(2.0))]),
            ),
        ]));
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        let layout = SystemLayout::from_phase_map(&HashMap::from([
            (Some("gas".to_string()), gas),
            (Some("liquid".to_string()), liquid),
        ]));
        let matrix = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]);

        let equations =
            build_symbolic_lagrange_equations(&layout, system.dG_sym_view(), &matrix, 300.0)
                .unwrap();
        assert_eq!(equations.len(), 2);
        assert_ne!(format!("{}", equations[0]), format!("{}", equations[1]));
        assert!(format!("{}", equations[0]).contains("Lambda0"));
        assert!(format!("{}", equations[1]).contains("Lambda1"));
    }

    #[test]
    fn symbolic_adapter_rejects_matrix_layout_mismatch() {
        let system = PhaseSystem::new_single_phase();
        let layout = SystemLayout::from_single_phase(&["A".to_string()]);
        let err = build_symbolic_lagrange_equations(
            &layout,
            system.dG_sym_view(),
            &DMatrix::zeros(2, 1),
            300.0,
        )
        .unwrap_err();
        assert!(err.to_string().contains("component rows"));
    }

    #[test]
    fn lagrange_adapters_reject_non_positive_or_non_finite_reference_temperature() {
        let layout = SystemLayout::from_single_phase(&["A".to_string()]);
        let matrix = DMatrix::from_row_slice(1, 1, &[1.0]);
        let system = PhaseSystem::new_single_phase();

        for reference_temperature in [0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let symbolic_error = build_symbolic_lagrange_equations(
                &layout,
                system.dG_sym_view(),
                &matrix,
                reference_temperature,
            )
            .unwrap_err();
            assert!(symbolic_error.to_string().contains("reference temperature"));

            let numeric_error = match build_numeric_lagrange_equations(
                &layout,
                HashMap::new(),
                matrix.clone(),
                reference_temperature,
            ) {
                Ok(_) => panic!("invalid reference temperature must be rejected"),
                Err(error) => error,
            };
            assert!(numeric_error.to_string().contains("reference temperature"));
        }
    }

    #[test]
    fn numeric_adapter_uses_global_component_indices_across_phases() {
        let mut gas = SubsData::new();
        gas.substances = vec!["H2O".to_string()];
        let mut liquid = SubsData::new();
        liquid.substances = vec!["H2O".to_string()];
        let layout = SystemLayout::from_phase_map(&HashMap::from([
            (Some("gas".to_string()), gas),
            (Some("liquid".to_string()), liquid),
        ]));
        let gibbs = HashMap::from([
            (
                Some("gas".to_string()),
                HashMap::from([(
                    "H2O".to_string(),
                    Box::new(|_, _, _| 1.0) as LegacyGibbsFunction,
                )]),
            ),
            (
                Some("liquid".to_string()),
                HashMap::from([(
                    "H2O".to_string(),
                    Box::new(|_, _, _| 2.0) as LegacyGibbsFunction,
                )]),
            ),
        ]);
        let equations = build_numeric_lagrange_equations(
            &layout,
            gibbs,
            DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 1.0]),
            300.0,
        )
        .unwrap();
        let result = equations(300.0, Some(vec![0.5, 0.5]), Some(1.0), vec![3.0, 4.0]);

        assert_eq!(result.len(), 2);
        assert!((result[0] - (3.0 + 1.0 / (R * 300.0))).abs() < 1e-12);
        assert!((result[1] - (4.0 + 2.0 / (R * 300.0))).abs() < 1e-12);
    }

    #[test]
    fn numeric_adapter_rejects_missing_gibbs_before_publishing_a_closure() {
        let layout = SystemLayout::from_single_phase(&["A".to_string()]);
        let error = match build_numeric_lagrange_equations(
            &layout,
            HashMap::new(),
            DMatrix::from_row_slice(1, 1, &[1.0]),
            300.0,
        ) {
            Ok(_) => panic!("missing Gibbs data must fail during closure construction"),
            Err(error) => error,
        };

        assert!(error.to_string().contains("phase Gibbs functions"));
    }

    #[test]
    fn numeric_adapter_runtime_guards_reject_malformed_solver_vectors() {
        let layout = SystemLayout::from_single_phase(&["A".to_string()]);
        let gibbs = HashMap::from([(
            None,
            HashMap::from([(
                "A".to_string(),
                Box::new(|_, _, _| 1.0) as LegacyGibbsFunction,
            )]),
        )]);
        let equations = build_numeric_lagrange_equations(
            &layout,
            gibbs,
            DMatrix::from_row_slice(1, 1, &[1.0]),
            300.0,
        )
        .expect("valid data must construct a closure");

        // The legacy closure cannot return `Result`; an empty vector is the
        // established non-panicking failure signal for malformed solver input.
        assert!(equations(300.0, Some(vec![]), Some(1.0), vec![1.0]).is_empty());
        assert!(equations(300.0, Some(vec![1.0]), Some(1.0), vec![]).is_empty());
        assert!(equations(f64::NAN, Some(vec![1.0]), Some(1.0), vec![1.0]).is_empty());
    }
}
