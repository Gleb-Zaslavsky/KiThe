//! Staged construction of symbolic thermodynamic properties and closures.
//!
//! These helpers prepare cloned phase payloads and return derived values; they
//! never publish cache state. `PhaseSystem` owns the commit, revision bump, and
//! cache-family replacement after a complete stage succeeds.

use std::collections::HashMap;
use std::sync::Arc;

use RustedSciThe::symbolic::symbolic_engine::Expr;

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::SubsDataResult;

use super::{PhaseFunction, PhaseSymbolicLayout, stage_phase_data_operation};

type StagedPhaseValues<T> = (
    HashMap<Option<String>, SubsData>,
    HashMap<Option<String>, HashMap<String, T>>,
);

/// Stages symbolic Gibbs expressions against the current mole-variable layout.
pub(crate) fn stage_symbolic_gibbs(
    phase_data: &HashMap<Option<String>, SubsData>,
    symbolic_layout: &PhaseSymbolicLayout,
    temperature: f64,
) -> SubsDataResult<StagedPhaseValues<Expr>> {
    let symbolic_vars = symbolic_layout.phase_variables().clone();
    stage_phase_data_operation(phase_data, |phase_name, subs_data| {
        let (phase_total, component_amounts) = symbolic_vars
            .get(phase_name)
            .cloned()
            .unwrap_or((None, None));
        subs_data.calculate_Gibbs_sym_one_phase(temperature, component_amounts, phase_total)
    })
}

/// Stages symbolic entropy expressions against the current mole-variable layout.
pub(crate) fn stage_symbolic_entropy(
    phase_data: &HashMap<Option<String>, SubsData>,
    symbolic_layout: &PhaseSymbolicLayout,
    temperature: f64,
) -> SubsDataResult<StagedPhaseValues<Expr>> {
    let symbolic_vars = symbolic_layout.phase_variables().clone();
    stage_phase_data_operation(phase_data, |phase_name, subs_data| {
        let (phase_total, component_amounts) = symbolic_vars
            .get(phase_name)
            .cloned()
            .unwrap_or((None, None));
        subs_data.calculate_S_sym_for_one_phase(temperature, component_amounts, phase_total)
    })
}

/// Stages Gibbs closures and converts legacy boxed functions into shareable arcs.
pub(crate) fn stage_gibbs_functions(
    phase_data: &HashMap<Option<String>, SubsData>,
    temperature: f64,
    pressure: f64,
) -> SubsDataResult<StagedPhaseValues<PhaseFunction>> {
    let (staged_data, values) = stage_phase_data_operation(phase_data, |_, subs_data| {
        subs_data.calculate_Gibbs_fun_one_phase(pressure, temperature)
    })?;
    Ok((staged_data, arcify_function_maps(values)))
}

/// Stages entropy closures and converts legacy boxed functions into shareable arcs.
pub(crate) fn stage_entropy_functions(
    phase_data: &HashMap<Option<String>, SubsData>,
    temperature: f64,
    pressure: f64,
) -> SubsDataResult<StagedPhaseValues<PhaseFunction>> {
    let (staged_data, values) = stage_phase_data_operation(phase_data, |_, subs_data| {
        subs_data.calculate_S_fun_for_one_phase(pressure, temperature)
    })?;
    Ok((staged_data, arcify_function_maps(values)))
}

fn arcify_function_maps(
    values: HashMap<
        Option<String>,
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>,
    >,
) -> HashMap<Option<String>, HashMap<String, PhaseFunction>> {
    values
        .into_iter()
        .map(|(phase, phase_values)| {
            let phase_values = phase_values
                .into_iter()
                .map(|(substance, function)| (substance, Arc::from(function)))
                .collect();
            (phase, phase_values)
        })
        .collect()
}
