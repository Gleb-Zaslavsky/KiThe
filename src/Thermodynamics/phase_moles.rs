//! Canonical mole-variable construction and mole-input normalization.
//!
//! This module translates legacy phase-keyed maps into snapshots aligned with
//! `SystemLayout`. It owns neither thermodynamic payloads nor cache lifecycle;
//! `PhaseSystem` supplies the layout and publishes the derived symbolic state.

use std::collections::HashMap;

use RustedSciThe::symbolic::symbolic_engine::Expr;

use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;

use super::{MoleNumberSnapshot, OrderedPhaseMoles, PhaseMoleNumbers, PhaseSymbolicLayout};

/// Builds canonical mole variables from the supplied phase/component layout.
pub(crate) fn build_indexed_mole_variables(
    layout: &SystemLayout,
) -> SubsDataResult<PhaseSymbolicLayout> {
    let mut all_component_vars = Vec::with_capacity(layout.component_count());
    let phase_count = layout.phases().len();
    let mut phase_totals = Vec::with_capacity(phase_count);
    let mut indexed_vars = HashMap::with_capacity(phase_count);
    let mut variables_by_substance = HashMap::with_capacity(phase_count);

    for phase in layout.phases() {
        let phase_name = phase.as_option().clone();
        let phase_index = layout
            .phase_index(phase)
            .ok_or_else(|| SubsDataError::MissingData {
                field: "phase index".to_string(),
                substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
            })?;
        let total = Expr::Var(format!("Np{phase_index}"));
        let components =
            layout
                .components_for_phase(phase)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase component layout".to_string(),
                    substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                })?;
        let mut phase_vars = Vec::with_capacity(components.len());
        let mut phase_substance_vars = HashMap::with_capacity(components.len());
        for (component_index, component) in components.iter().enumerate() {
            let amount = Expr::Var(format!("n{phase_index}_{component_index}"));
            all_component_vars.push(amount.clone());
            phase_vars.push(amount.clone());
            phase_substance_vars.insert(component.substance.clone(), amount);
        }
        phase_totals.push(total.clone());
        indexed_vars.insert(phase_name.clone(), (Some(total), Some(phase_vars)));
        variables_by_substance.insert(phase_name, phase_substance_vars);
    }

    Ok(PhaseSymbolicLayout::new(
        indexed_vars,
        all_component_vars,
        phase_totals,
        variables_by_substance,
    ))
}

/// Normalizes sparse phase-keyed mole input into the supplied canonical layout.
pub(crate) fn normalize_mole_numbers(
    layout: SystemLayout,
    non_zero_number_of_moles: HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
) -> SubsDataResult<MoleNumberSnapshot> {
    let mut phase_amounts = non_zero_number_of_moles
        .into_iter()
        .map(|(phase, (total_amount, component_amounts))| {
            (
                phase,
                PhaseMoleNumbers::new(total_amount, component_amounts),
            )
        })
        .collect::<HashMap<_, _>>();
    let mut ordered_phase_amounts = HashMap::with_capacity(layout.phases().len());
    let mut totals_by_substance = HashMap::new();

    for phase in layout.phases() {
        let phase_name = phase.as_option().clone();
        let Some(phase_amounts) = phase_amounts.get_mut(&phase_name) else {
            continue;
        };
        let Some(amounts) = phase_amounts.component_amounts.as_mut() else {
            continue;
        };
        let components =
            layout
                .components_for_phase(phase)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase component layout".to_string(),
                    substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                })?;
        for component in components {
            amounts.entry(component.substance.clone()).or_insert(0.0);
        }
        let ordered_amounts = components
            .iter()
            .map(|component| *amounts.get(&component.substance).unwrap_or(&0.0))
            .collect::<Vec<_>>();
        for (substance, amount) in amounts.iter() {
            *totals_by_substance.entry(substance.clone()).or_insert(0.0) += amount;
        }
        ordered_phase_amounts.insert(
            phase_name,
            OrderedPhaseMoles {
                total_amount: phase_amounts.total_amount,
                component_amounts: ordered_amounts,
            },
        );
    }

    Ok(MoleNumberSnapshot {
        layout,
        phase_amounts,
        ordered_phase_amounts,
        total_amounts_by_substance: totals_by_substance,
    })
}
