//! Typed inputs and pure evaluation helpers for phase thermodynamics.
//!
//! MAIN DATA STRUCTURES:
//!
//! ThermoEvaluationConditions - shared temperature and pressure for a numeric evaluation
//! PhaseComposition - typed phase composition payload for a numeric evaluation
//!
//!
//!IERARCHY:
//!                              NumericThermoCacheContext
//!                             /
//!                         PhaseEvaluationRequest   
//!                         /                     \
//! ThermoEvaluationConditions                   PhaseComposition
//!
//! This module owns the request boundary between a resolved phase system and
//! raw `SubsData` calculators.  It intentionally does not own cache state or
//! facade types: the same validated request can be evaluated by either a
//! single-phase or multi-phase view.
//!

use std::collections::{HashMap, HashSet};

use crate::Thermodynamics::User_substances::{Phases, SubsData};
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::{SystemLayout, ordered_phase_keys};

/// Typed phase composition payload used internally instead of anonymous tuples.
/// `total_amount` is the phase total (`Np`), while `component_amounts` keeps
/// the ordered mole vector for the phase components when it is available.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseComposition {
    total_amount: Option<f64>,
    component_amounts: Option<Vec<f64>>,
}

impl PhaseComposition {
    pub fn new(total_amount: Option<f64>, component_amounts: Option<Vec<f64>>) -> Self {
        Self {
            total_amount,
            component_amounts,
        }
    }

    pub fn try_new_numeric(
        total_amount: Option<f64>,
        component_amounts: Option<Vec<f64>>,
        phase_name: &Option<String>,
        operation: &str,
    ) -> SubsDataResult<Self> {
        let composition = Self::new(total_amount, component_amounts);
        composition.validate_numeric(phase_name, operation)?;
        Ok(composition)
    }

    pub fn into_legacy(self) -> (Option<f64>, Option<Vec<f64>>) {
        (self.total_amount, self.component_amounts)
    }

    pub fn total_amount(&self) -> Option<f64> {
        self.total_amount
    }

    pub fn component_amounts(&self) -> Option<&[f64]> {
        self.component_amounts.as_deref()
    }

    pub(crate) fn component_amounts_owned(&self) -> Option<Vec<f64>> {
        self.component_amounts.clone()
    }

    /// Validates a numeric phase composition before Gibbs/entropy evaluation.
    /// The phase total and every component amount must be finite and
    /// non-negative, and the total must match the summed component amounts.
    pub fn validate_numeric(
        &self,
        phase_name: &Option<String>,
        operation: &str,
    ) -> SubsDataResult<()> {
        let phase_label = phase_name.clone().unwrap_or_else(|| "None".to_string());
        let total = self.total_amount.ok_or_else(|| {
            SubsDataError::calculation_failed(
                phase_label.clone(),
                operation,
                "missing phase total amount",
            )
        })?;
        if !total.is_finite() || total <= 0.0 {
            return Err(SubsDataError::calculation_failed(
                phase_label.clone(),
                operation,
                format!("invalid phase total amount {}", total),
            ));
        }

        let components = self.component_amounts.as_ref().ok_or_else(|| {
            SubsDataError::calculation_failed(
                phase_label.clone(),
                operation,
                "missing component amount vector",
            )
        })?;
        if components.is_empty() {
            return Err(SubsDataError::calculation_failed(
                phase_label.clone(),
                operation,
                "empty component amount vector",
            ));
        }

        let mut sum = 0.0;
        for value in components {
            if !value.is_finite() || *value < 0.0 {
                return Err(SubsDataError::calculation_failed(
                    phase_label.clone(),
                    operation,
                    format!("invalid component amount {}", value),
                ));
            }
            sum += value;
        }

        let tolerance = 1e-9_f64.max(total.abs() * 1e-9_f64);
        if (sum - total).abs() > tolerance {
            return Err(SubsDataError::calculation_failed(
                phase_label,
                operation,
                format!("phase total {} does not match component sum {}", total, sum),
            ));
        }
        Ok(())
    }
}

/// Physical conditions shared by a numeric phase-property evaluation.
///
/// Keeping temperature and pressure together prevents validating one input
/// while silently forwarding an invalid companion value to a calculator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoEvaluationConditions {
    temperature: f64,
    pressure: f64,
}

impl ThermoEvaluationConditions {
    pub fn new(temperature: f64, pressure: f64) -> SubsDataResult<Self> {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(SubsDataError::InvalidTemperature(temperature));
        }
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(SubsDataError::InvalidPhysicalValue {
                field: "pressure".to_string(),
                value: pressure,
            });
        }
        Ok(Self {
            temperature,
            pressure,
        })
    }

    pub fn temperature(self) -> f64 {
        self.temperature
    }

    pub fn pressure(self) -> f64 {
        self.pressure
    }
}

/// Complete input to one numeric phase-property evaluation.
///
/// Temperature, pressure, and phase compositions are one logical request:
/// separating them made it too easy to evaluate with one composition and
/// publish cache provenance for another.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEvaluationRequest {
    conditions: ThermoEvaluationConditions,
    compositions: HashMap<Option<String>, PhaseComposition>,
}

impl PhaseEvaluationRequest {
    pub fn new(
        conditions: ThermoEvaluationConditions,
        compositions: HashMap<Option<String>, PhaseComposition>,
    ) -> Self {
        Self {
            conditions,
            compositions,
        }
    }

    pub fn from_numeric(
        temperature: f64,
        pressure: f64,
        compositions: HashMap<Option<String>, PhaseComposition>,
    ) -> SubsDataResult<Self> {
        Ok(Self::new(
            ThermoEvaluationConditions::new(temperature, pressure)?,
            compositions,
        ))
    }

    pub fn conditions(&self) -> ThermoEvaluationConditions {
        self.conditions
    }

    pub fn compositions(&self) -> &HashMap<Option<String>, PhaseComposition> {
        &self.compositions
    }
}

/// Provenance for a cached numeric phase-property result.
///
/// A numeric result is reusable only for the exact physical conditions and
/// phase compositions that produced it, and only while the resolved phase
/// configuration is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct NumericThermoCacheContext {
    request: PhaseEvaluationRequest,
    configuration_revision: usize,
}

impl NumericThermoCacheContext {
    pub(crate) fn new(request: PhaseEvaluationRequest, configuration_revision: usize) -> Self {
        Self {
            request,
            configuration_revision,
        }
    }

    pub fn conditions(&self) -> ThermoEvaluationConditions {
        self.request.conditions()
    }

    pub fn compositions(&self) -> &HashMap<Option<String>, PhaseComposition> {
        self.request.compositions()
    }

    pub fn request(&self) -> &PhaseEvaluationRequest {
        &self.request
    }

    pub fn configuration_revision(&self) -> usize {
        self.configuration_revision
    }

    pub(crate) fn matches_request(
        &self,
        request: &PhaseEvaluationRequest,
        configuration_revision: usize,
    ) -> bool {
        self.request == *request && self.configuration_revision == configuration_revision
    }

    pub(crate) fn matches_parts(
        &self,
        conditions: ThermoEvaluationConditions,
        compositions: &HashMap<Option<String>, PhaseComposition>,
        configuration_revision: usize,
    ) -> bool {
        self.request.conditions() == conditions
            && self.request.compositions() == compositions
            && self.configuration_revision == configuration_revision
    }
}

pub(crate) fn typed_phase_composition_for<'a>(
    compositions: &'a HashMap<Option<String>, PhaseComposition>,
    phase_name: &Option<String>,
) -> SubsDataResult<&'a PhaseComposition> {
    compositions
        .get(phase_name)
        .ok_or_else(|| SubsDataError::MissingData {
            field: "phase data".to_string(),
            substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
        })
}

/// Validates the complete numeric request before any phase calculator runs.
///
/// The supplied phase keys must exactly match the resolved system and each
/// amount vector must follow that phase's declared component order. This
/// front-loads failure so a later malformed phase cannot alter an earlier one.
pub(crate) fn validate_phase_evaluation_request(
    expected_components: impl IntoIterator<Item = (Option<String>, usize)>,
    conditions: ThermoEvaluationConditions,
    compositions: &HashMap<Option<String>, PhaseComposition>,
    operation: &str,
) -> SubsDataResult<()> {
    let expected = expected_components.into_iter().collect::<HashMap<_, _>>();
    let supplied = compositions.keys().cloned().collect::<HashSet<_>>();
    let declared = expected.keys().cloned().collect::<HashSet<_>>();
    if supplied != declared {
        return Err(SubsDataError::calculation_failed(
            "phase system",
            operation,
            "composition phases must exactly match the resolved phase set",
        ));
    }

    for (phase, component_count) in expected {
        let composition = typed_phase_composition_for(compositions, &phase)?;
        composition.validate_numeric(&phase, operation)?;
        let actual_count = composition.component_amounts().map_or(0, <[f64]>::len);
        if actual_count != component_count {
            return Err(SubsDataError::calculation_failed(
                phase.clone().unwrap_or_else(|| "None".to_string()),
                operation,
                format!(
                    "component amount vector has length {}; expected {}",
                    actual_count, component_count
                ),
            ));
        }
    }

    let _ = conditions;
    Ok(())
}

/// Read-only bridge over the still transitional raw `SubsData` storage.
///
/// `OnePhase` keeps a direct payload while `PhaseOrSolution` currently keeps
/// a phase-keyed map because equilibrium consumers still mutate those legacy
/// fields. New algorithms receive this view and never publish calculator
/// state back into either facade.
#[derive(Clone, Copy)]
pub(crate) enum PhaseDataView<'a> {
    Single(&'a SubsData),
    Multi(&'a HashMap<Option<String>, SubsData>),
}

impl<'a> PhaseDataView<'a> {
    pub(crate) fn single(data: &'a SubsData) -> Self {
        Self::Single(data)
    }

    pub(crate) fn multi(data: &'a HashMap<Option<String>, SubsData>) -> Self {
        Self::Multi(data)
    }

    pub(crate) fn layout(self) -> SystemLayout {
        match self {
            Self::Single(data) => SystemLayout::from_single_phase(&data.substances),
            Self::Multi(data) => SystemLayout::from_phase_map(data),
        }
    }

    fn expected_component_counts(self) -> Vec<(Option<String>, usize)> {
        match self {
            Self::Single(data) => vec![(None, data.substances.len())],
            Self::Multi(data) => data
                .iter()
                .map(|(phase, data)| (phase.clone(), data.substances.len()))
                .collect(),
        }
    }

    /// Evaluation receives a private mutable copy because `SubsData` still
    /// prepares coefficient caches internally. This preserves the public
    /// `&self` evaluation contract during the migration.
    fn to_owned_phase_map(self) -> HashMap<Option<String>, SubsData> {
        match self {
            Self::Single(data) => HashMap::from([(None, data.clone())]),
            Self::Multi(data) => data.clone(),
        }
    }

    /// Rejects zero gas-component amounts before the raw calculators evaluate
    /// their chemical-potential correction `ln(n_i / Np)`.
    ///
    /// A zero amount remains meaningful for a condensed component, but a
    /// component-wise ideal-gas chemical potential is undefined at zero mole
    /// fraction.  The legacy calculator also treats a missing phase entry as
    /// gas, so this validation deliberately preserves that established rule.
    fn validate_gas_composition_domain(
        self,
        compositions: &HashMap<Option<String>, PhaseComposition>,
        operation: &str,
    ) -> SubsDataResult<()> {
        let validate_phase = |phase_name: &Option<String>, subs_data: &SubsData| {
            let composition = typed_phase_composition_for(compositions, phase_name)?;
            let component_amounts = composition.component_amounts().ok_or_else(|| {
                SubsDataError::calculation_failed(
                    phase_name.clone().unwrap_or_else(|| "None".to_string()),
                    operation,
                    "missing component amount vector",
                )
            })?;

            for (substance, amount) in subs_data.substances.iter().zip(component_amounts) {
                let uses_ideal_gas_correction = matches!(
                    subs_data.map_of_phases.get(substance),
                    None | Some(Some(Phases::Gas))
                );
                if uses_ideal_gas_correction && *amount == 0.0 {
                    return Err(SubsDataError::calculation_failed(
                        substance.clone(),
                        operation,
                        format!(
                            "zero component amount in phase {:?} is outside the ideal-gas chemical-potential domain",
                            phase_name
                        ),
                    ));
                }
            }
            Ok(())
        };

        match self {
            Self::Single(data) => validate_phase(&None, data),
            Self::Multi(data) => {
                for phase_name in ordered_phase_keys(data) {
                    let subs_data =
                        data.get(&phase_name)
                            .ok_or_else(|| SubsDataError::MissingData {
                                field: "resolved phase data".to_string(),
                                substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                            })?;
                    validate_phase(&phase_name, subs_data)?;
                }
                Ok(())
            }
        }
    }
}

/// Numeric property family evaluated by the shared raw-data bridge.
#[derive(Clone, Copy)]
pub(crate) enum NumericPhaseProperty {
    Gibbs,
    Entropy,
}

impl NumericPhaseProperty {
    fn operation_name(self) -> &'static str {
        match self {
            Self::Gibbs => "Gibbs free energy",
            Self::Entropy => "entropy",
        }
    }
}

/// Evaluates one numeric property over either a single or multi-phase view.
/// Validation intentionally happens before the first calculator runs.
pub(crate) fn evaluate_numeric_phase_property(
    data: PhaseDataView<'_>,
    request: &PhaseEvaluationRequest,
    property: NumericPhaseProperty,
) -> SubsDataResult<HashMap<Option<String>, HashMap<String, f64>>> {
    validate_phase_evaluation_request(
        data.expected_component_counts(),
        request.conditions(),
        request.compositions(),
        property.operation_name(),
    )?;
    data.validate_gas_composition_domain(request.compositions(), property.operation_name())?;

    let mut phase_data = data.to_owned_phase_map();
    let mut values = HashMap::with_capacity(phase_data.len());
    for phase_name in ordered_phase_keys(&phase_data) {
        let subsdata =
            phase_data
                .get_mut(&phase_name)
                .ok_or_else(|| SubsDataError::MissingData {
                    field: "resolved phase data".to_string(),
                    substance: phase_name.clone().unwrap_or_else(|| "None".to_string()),
                })?;
        let composition = typed_phase_composition_for(request.compositions(), &phase_name)?;
        let phase_values = match property {
            NumericPhaseProperty::Gibbs => subsdata.calc_dG_for_one_phase(
                request.conditions().pressure(),
                request.conditions().temperature(),
                composition.component_amounts_owned(),
                composition.total_amount(),
            )?,
            NumericPhaseProperty::Entropy => subsdata.calculate_S_for_one_phase(
                request.conditions().pressure(),
                request.conditions().temperature(),
                composition.component_amounts_owned(),
                composition.total_amount(),
            )?,
        };
        values.insert(phase_name, phase_values);
    }
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request_with_zero_first_component() -> PhaseEvaluationRequest {
        PhaseEvaluationRequest::new(
            ThermoEvaluationConditions::new(300.0, 101_325.0).unwrap(),
            HashMap::from([(
                None,
                PhaseComposition::try_new_numeric(Some(1.0), Some(vec![0.0, 1.0]), &None, "test")
                    .unwrap(),
            )]),
        )
    }

    fn phase_data(phase: Option<Phases>) -> SubsData {
        let mut data = SubsData::new();
        data.substances = vec!["A".to_string(), "B".to_string()];
        data.map_of_phases.insert("A".to_string(), phase);
        data.map_of_phases.insert("B".to_string(), phase);
        data
    }

    #[test]
    fn zero_gas_component_is_rejected_before_numeric_calculation() {
        let data = phase_data(Some(Phases::Gas));
        let error = evaluate_numeric_phase_property(
            PhaseDataView::single(&data),
            &request_with_zero_first_component(),
            NumericPhaseProperty::Gibbs,
        )
        .expect_err("zero ideal-gas amount must not reach ln(n_i / Np)");

        assert!(
            error
                .to_string()
                .contains("ideal-gas chemical-potential domain"),
            "unexpected error: {error}"
        );
    }

    #[test]
    fn zero_condensed_component_remains_valid_at_the_phase_boundary() {
        let data = phase_data(Some(Phases::Solid));
        PhaseDataView::single(&data)
            .validate_gas_composition_domain(
                request_with_zero_first_component().compositions(),
                "test",
            )
            .expect("zero condensed amount is outside the ideal-gas domain and remains valid");
    }
}
