//! Pure element-matrix assembly from resolved phase payloads.
//!
//! The returned rows follow `SystemLayout` component order. The helper clones
//! each `SubsData` payload because lower-level composition extraction may
//! prepare transient calculator state; read-only phase queries must not change
//! the engine's stored payloads.

use std::collections::{HashMap, HashSet};

use nalgebra::DMatrix;

use crate::Thermodynamics::User_substances::SubsData;
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;

/// Returns a component-by-element matrix, molar masses, and sorted element labels.
pub(crate) fn element_composition_and_molar_mass(
    phase_data: &HashMap<Option<String>, SubsData>,
    layout: &SystemLayout,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> SubsDataResult<(DMatrix<f64>, HashMap<String, f64>, Vec<String>)> {
    let mut component_compositions = Vec::with_capacity(layout.component_count());
    let mut elements = HashSet::new();
    let mut molar_masses = HashMap::new();

    for phase in layout.phases() {
        let phase_name = phase.as_option();
        let phase_data = phase_data
            .get(phase_name)
            .ok_or_else(|| SubsDataError::MissingData {
                field: "phase data".to_string(),
                substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
            })?;
        let mut local_phase_data = phase_data.clone();
        let (phase_molar_masses, compositions, phase_elements) =
            SubsData::calculate_elem_composition_and_molar_mass_local(
                &mut local_phase_data,
                groups.clone(),
            )
            .map_err(|error| SubsDataError::MatrixOperationFailed(error.to_string()))?;

        let phase_components =
            layout
                .components_for_phase(phase)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "phase component layout".to_string(),
                    substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                })?;
        for component_index in 0..phase_components.len() {
            component_compositions.push(compositions.get(component_index).cloned().ok_or_else(
                || {
                    SubsDataError::MatrixOperationFailed(
                        "component composition index out of bounds".to_string(),
                    )
                },
            )?);
        }
        elements.extend(phase_elements);
        molar_masses.extend(phase_molar_masses);
    }

    let mut unique_elements = elements.into_iter().collect::<Vec<_>>();
    unique_elements.sort();
    let mut matrix = DMatrix::zeros(unique_elements.len(), component_compositions.len());
    for (component_index, composition) in component_compositions.iter().enumerate() {
        for (element_index, element) in unique_elements.iter().enumerate() {
            if let Some(count) = composition.get(element) {
                matrix[(element_index, component_index)] = *count;
            }
        }
    }

    Ok((matrix.transpose(), molar_masses, unique_elements))
}
