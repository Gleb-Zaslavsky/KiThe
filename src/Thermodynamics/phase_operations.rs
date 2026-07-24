//! Transactional operations over cloned phase payloads.
//!
//! Phase-local `SubsData` preparation can fail. This helper guarantees that a
//! caller receives either every updated phase payload and every result, or an
//! error with the original map untouched. `PhaseSystem` alone decides when a
//! successful staged map becomes visible and which caches must be invalidated.

use std::collections::HashMap;

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::SubsDataResult;

/// Runs a fallible operation over cloned phase payloads as one transaction.
pub(crate) fn stage_phase_data_operation<R>(
    phase_data: &HashMap<Option<String>, SubsData>,
    mut update: impl FnMut(&Option<String>, &mut SubsData) -> SubsDataResult<R>,
) -> SubsDataResult<(
    HashMap<Option<String>, SubsData>,
    HashMap<Option<String>, R>,
)> {
    let mut updated_phase_data = phase_data.clone();
    let mut results = HashMap::with_capacity(updated_phase_data.len());
    for (phase_name, subs_data) in updated_phase_data.iter_mut() {
        results.insert(phase_name.clone(), update(phase_name, subs_data)?);
    }
    Ok((updated_phase_data, results))
}
