//! Typed presets for common thermodynamics workflows.
//!
//! The goal is to keep scenario-specific setup out of `User_substances.rs`
//! while still reusing its existing setters and lookup policy machinery.

use crate::Thermodynamics::prelude::{
    DataType, LibraryId, LibraryPriority, PropertyKind, SubsData, SubsDataResult,
};
use std::collections::HashMap;

/// High-level thermodynamics workflow family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermoScenarioKind {
    /// Stirred tank workflows usually need thermochemistry only.
    StirredTank,
    /// Plug-flow gas workflows often need thermochemistry and transport.
    PlugFlowGas,
    /// Plug-flow condensed workflows need thermochemistry and conductivity only.
    PlugFlowCond,
}

/// Typed lookup policy for a thermodynamics preset.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThermoLookupPolicy {
    pub priority_libraries: Vec<LibraryId>,
    pub permitted_libraries: Vec<LibraryId>,
    pub explicit_instructions: Vec<(String, LibraryId)>,
    pub temperature_interval: Option<(f64, f64)>,
}

impl ThermoLookupPolicy {
    /// Creates an empty lookup policy.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a priority library to the lookup policy.
    pub fn with_priority_library(mut self, library: LibraryId) -> Self {
        self.priority_libraries.push(library);
        self
    }

    /// Adds a permitted fallback library to the lookup policy.
    pub fn with_permitted_library(mut self, library: LibraryId) -> Self {
        self.permitted_libraries.push(library);
        self
    }

    /// Adds an explicit substance-to-library mapping to the lookup policy.
    pub fn with_explicit_instruction(mut self, substance: String, library: LibraryId) -> Self {
        self.explicit_instructions.push((substance, library));
        self
    }

    /// Sets the temperature interval carried by the lookup policy.
    pub fn with_temperature_interval(mut self, lower: f64, upper: f64) -> Self {
        self.temperature_interval = Some((lower, upper));
        self
    }

    /// Returns `true` if the policy does not impose any lookup constraints.
    pub fn is_empty(&self) -> bool {
        self.priority_libraries.is_empty()
            && self.permitted_libraries.is_empty()
            && self.explicit_instructions.is_empty()
            && self.temperature_interval.is_none()
    }
}

/// Output-style policy for scenario previews and downstream integration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermoOutputPolicy {
    /// Return only the requested thermodynamic properties.
    PropertiesOnly,
    /// Return thermodynamic properties and transport data.
    PropertiesAndTransport,
    /// Return thermodynamic properties plus a preview report.
    PropertiesAndPreviewReport,
    /// Return thermodynamic properties plus a write-workflow report.
    PropertiesAndWriteWorkflow,
}

impl ThermoOutputPolicy {
    /// Returns a short human-readable label.
    pub const fn label(self) -> &'static str {
        match self {
            Self::PropertiesOnly => "properties_only",
            Self::PropertiesAndTransport => "properties_and_transport",
            Self::PropertiesAndPreviewReport => "properties_and_preview_report",
            Self::PropertiesAndWriteWorkflow => "properties_and_write_workflow",
        }
    }
}

/// Typed scenario preset for common thermodynamics workflows.
///
/// This is intentionally a configuration object rather than a solver wrapper:
/// it describes what the caller wants, and `apply_to` maps that request onto
/// the existing `SubsData` setters.
#[derive(Debug, Clone)]
pub struct ThermoScenarioPreset {
    pub kind: ThermoScenarioKind,
    pub substances: Vec<String>,
    pub lookup_policy: ThermoLookupPolicy,
    pub output_policy: ThermoOutputPolicy,
    pub temperature: Option<f64>,
    pub pressure: Option<f64>,
    pub molar_masses: Option<HashMap<String, f64>>,
    pub needs_molar_masses: bool,
    pub needs_thermochemistry: bool,
    pub needs_lambda: bool,
    pub needs_viscosity: bool,
    pub needs_diffusion: bool,
}

impl ThermoScenarioPreset {
    fn new(kind: ThermoScenarioKind, substances: Vec<String>) -> Self {
        Self {
            kind,
            substances,
            lookup_policy: ThermoLookupPolicy::new(),
            output_policy: ThermoOutputPolicy::PropertiesOnly,
            temperature: None,
            pressure: None,
            molar_masses: None,
            needs_molar_masses: true,
            needs_thermochemistry: true,
            needs_lambda: false,
            needs_viscosity: false,
            needs_diffusion: false,
        }
    }

    /// Preset for stirred-tank calculations.
    pub fn stirred_tank(substances: Vec<String>) -> Self {
        Self::new(ThermoScenarioKind::StirredTank, substances)
    }

    /// Preset for plug-flow gas calculations.
    pub fn plug_flow_gas(substances: Vec<String>) -> Self {
        let mut preset = Self::new(ThermoScenarioKind::PlugFlowGas, substances);
        preset.output_policy = ThermoOutputPolicy::PropertiesAndTransport;
        preset.needs_lambda = true;
        preset.needs_viscosity = true;
        preset.needs_diffusion = true;
        preset
    }

    /// Preset for plug-flow condensed calculations.
    pub fn plug_flow_cond(substances: Vec<String>) -> Self {
        let mut preset = Self::new(ThermoScenarioKind::PlugFlowCond, substances);
        preset.output_policy = ThermoOutputPolicy::PropertiesAndTransport;
        preset.needs_lambda = true;
        preset.needs_viscosity = false;
        preset.needs_diffusion = false;
        preset
    }

    /// Preset for plug-flow gas calculations without diffusion.
    pub fn plug_flow_gas_without_diffusion(substances: Vec<String>) -> Self {
        let mut preset = Self::plug_flow_gas(substances);
        preset.needs_diffusion = false;
        preset
    }

    /// Sets a representative temperature for the preset.
    pub fn with_temperature(mut self, temperature: f64) -> Self {
        self.temperature = Some(temperature);
        self
    }

    /// Sets the temperature interval used by the scenario description.
    pub fn with_temperature_interval(mut self, lower: f64, upper: f64) -> Self {
        self.lookup_policy = self.lookup_policy.with_temperature_interval(lower, upper);
        self
    }

    /// Sets the pressure for the preset.
    pub fn with_pressure(mut self, pressure: f64) -> Self {
        self.pressure = Some(pressure);
        self
    }

    /// Sets the molar-mass map for the preset.
    pub fn with_molar_masses(mut self, molar_masses: HashMap<String, f64>) -> Self {
        self.molar_masses = Some(molar_masses);
        self
    }

    /// Adds a priority library to the scenario.
    pub fn with_priority_library(mut self, library: LibraryId) -> Self {
        self.lookup_policy = self.lookup_policy.with_priority_library(library);
        self
    }

    /// Adds a permitted fallback library to the scenario.
    pub fn with_permitted_library(mut self, library: LibraryId) -> Self {
        self.lookup_policy = self.lookup_policy.with_permitted_library(library);
        self
    }

    /// Adds an explicit substance-to-library mapping to the scenario.
    pub fn with_explicit_instruction(mut self, substance: String, library: LibraryId) -> Self {
        self.lookup_policy = self
            .lookup_policy
            .with_explicit_instruction(substance, library);
        self
    }

    /// Sets the output policy for the preset.
    pub fn with_output_policy(mut self, output_policy: ThermoOutputPolicy) -> Self {
        self.output_policy = output_policy;
        self
    }

    /// Returns the canonical property set requested by the scenario.
    pub fn requested_properties(&self) -> Vec<PropertyKind> {
        let mut properties = vec![PropertyKind::Cp, PropertyKind::dH, PropertyKind::dS];
        if self.needs_lambda {
            properties.push(PropertyKind::Lambda);
        }
        if self.needs_viscosity {
            properties.push(PropertyKind::Visc);
        }
        properties
    }

    /// Returns the representation-specific data keys that the caller will
    /// typically need for this scenario.
    pub fn requested_data_types(&self) -> Vec<DataType> {
        self.requested_data_request().numeric_data_types
    }

    /// Returns the canonical data request needed by solver-oriented consumers.
    pub fn requested_data_request(&self) -> ThermoRequestedData {
        let properties = self.requested_properties();
        let numeric_data_types = properties
            .iter()
            .copied()
            .map(PropertyKind::numeric_type)
            .collect::<Vec<_>>();
        let function_data_types = properties
            .iter()
            .copied()
            .map(PropertyKind::function_type)
            .collect::<Vec<_>>();
        let symbolic_data_types = properties
            .iter()
            .copied()
            .map(PropertyKind::symbolic_type)
            .collect::<Vec<_>>();

        ThermoRequestedData {
            substances: self.substances.clone(),
            property_kinds: properties,
            numeric_data_types,
            function_data_types,
            symbolic_data_types,
            requires_molar_masses: self.needs_molar_masses,
        }
    }

    /// Returns the preset as human-readable request lines.
    pub fn requested_inputs(&self) -> Vec<String> {
        self.summary_rows()
            .into_iter()
            .map(|(key, value)| format!("{key}: {value}"))
            .collect()
    }

    /// Returns a compact summary of the preset for GUI preview tables.
    pub fn summary_rows(&self) -> Vec<(String, String)> {
        let mut rows = vec![
            ("scenario".to_string(), format!("{:?}", self.kind)),
            (
                "output_policy".to_string(),
                self.output_policy.label().to_string(),
            ),
            ("substances".to_string(), self.substances.join(", ")),
            (
                "requested_properties".to_string(),
                self.requested_properties()
                    .into_iter()
                    .map(|property| format!("{:?}", property))
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            (
                "requested_data_types".to_string(),
                self.requested_data_types()
                    .into_iter()
                    .map(|data_type| format!("{:?}", data_type))
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            (
                "requested_function_data_types".to_string(),
                self.requested_data_request()
                    .function_data_types
                    .iter()
                    .map(|data_type| format!("{:?}", data_type))
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            (
                "requested_symbolic_data_types".to_string(),
                self.requested_data_request()
                    .symbolic_data_types
                    .iter()
                    .map(|data_type| format!("{:?}", data_type))
                    .collect::<Vec<_>>()
                    .join(", "),
            ),
            (
                "needs_thermochemistry".to_string(),
                self.needs_thermochemistry.to_string(),
            ),
            (
                "needs_molar_masses".to_string(),
                self.needs_molar_masses.to_string(),
            ),
            ("needs_lambda".to_string(), self.needs_lambda.to_string()),
            (
                "needs_viscosity".to_string(),
                self.needs_viscosity.to_string(),
            ),
            (
                "needs_diffusion".to_string(),
                self.needs_diffusion.to_string(),
            ),
        ];

        if let Some((lower, upper)) = self.lookup_policy.temperature_interval {
            rows.push((
                "temperature_interval".to_string(),
                format!("{lower}..{upper}"),
            ));
        }

        if let Some(temperature) = self.temperature {
            rows.push(("temperature".to_string(), temperature.to_string()));
        }

        if let Some(pressure) = self.pressure {
            rows.push(("pressure".to_string(), pressure.to_string()));
        }

        if let Some(molar_masses) = &self.molar_masses {
            rows.push(("molar_masses".to_string(), format!("{:?}", molar_masses)));
        }

        if !self.lookup_policy.priority_libraries.is_empty() {
            rows.push((
                "priority_libraries".to_string(),
                self.lookup_policy
                    .priority_libraries
                    .iter()
                    .map(|library| library.canonical_name())
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }

        if !self.lookup_policy.permitted_libraries.is_empty() {
            rows.push((
                "permitted_libraries".to_string(),
                self.lookup_policy
                    .permitted_libraries
                    .iter()
                    .map(|library| library.canonical_name())
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }

        if !self.lookup_policy.explicit_instructions.is_empty() {
            rows.push((
                "explicit_instructions".to_string(),
                self.lookup_policy
                    .explicit_instructions
                    .iter()
                    .map(|(substance, library)| {
                        format!("{substance}=>{}", library.canonical_name())
                    })
                    .collect::<Vec<_>>()
                    .join(", "),
            ));
        }

        rows
    }

    /// Returns a typed preview report that can be rendered by a GUI or CLI.
    pub fn preview_report(&self) -> ThermoPresetPreviewReport {
        ThermoPresetPreviewReport {
            kind: self.kind,
            output_policy: self.output_policy,
            rows: self.summary_rows(),
        }
    }

    /// Applies the preset to an existing `SubsData` instance.
    pub fn apply_to(&self, subs_data: &mut SubsData) -> SubsDataResult<()> {
        subs_data.set_substances(self.substances.clone());

        if !self.lookup_policy.priority_libraries.is_empty() {
            subs_data.set_multiple_library_priorities_ids(
                self.lookup_policy.priority_libraries.clone(),
                LibraryPriority::Priority,
            );
        }

        if !self.lookup_policy.permitted_libraries.is_empty() {
            subs_data.set_multiple_library_priorities_ids(
                self.lookup_policy.permitted_libraries.clone(),
                LibraryPriority::Permitted,
            );
        }

        for (substance, library) in &self.lookup_policy.explicit_instructions {
            subs_data.set_explicit_search_instruction(substance.clone(), *library);
        }

        if let Some((lower, upper)) = self.lookup_policy.temperature_interval {
            subs_data.set_T_range_for_all_thermo(lower, upper)?;
        }

        if let Some(temperature) = self.temperature {
            subs_data.set_T(temperature)?;
        }

        if let Some(pressure) = self.pressure {
            subs_data.set_P(pressure, None)?;
        }

        if let Some(molar_masses) = &self.molar_masses {
            subs_data.set_M(molar_masses.clone(), None)?;
        }

        Ok(())
    }
}

/// Preview report for a thermodynamics preset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermoPresetPreviewReport {
    pub kind: ThermoScenarioKind,
    pub output_policy: ThermoOutputPolicy,
    pub rows: Vec<(String, String)>,
}

impl ThermoPresetPreviewReport {
    /// Returns the preview rows in a printable form.
    pub fn lines(&self) -> Vec<String> {
        self.rows
            .iter()
            .map(|(key, value)| format!("{key}: {value}"))
            .collect()
    }
}

/// Structured data request emitted by a thermodynamics preset.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermoRequestedData {
    pub substances: Vec<String>,
    pub property_kinds: Vec<PropertyKind>,
    pub numeric_data_types: Vec<DataType>,
    pub function_data_types: Vec<DataType>,
    pub symbolic_data_types: Vec<DataType>,
    pub requires_molar_masses: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stirred_tank_preset_requests_thermochemistry_only() {
        let preset = ThermoScenarioPreset::stirred_tank(vec!["CO".to_string()]);
        assert_eq!(preset.kind, ThermoScenarioKind::StirredTank);
        assert!(preset.needs_thermochemistry);
        assert!(!preset.needs_lambda);
        assert!(!preset.needs_viscosity);
        assert!(!preset.needs_diffusion);
        assert_eq!(preset.output_policy, ThermoOutputPolicy::PropertiesOnly);
        assert_eq!(
            preset.requested_properties(),
            vec![PropertyKind::Cp, PropertyKind::dH, PropertyKind::dS]
        );
    }

    #[test]
    fn plug_flow_gas_preset_requests_transport() {
        let preset = ThermoScenarioPreset::plug_flow_gas(vec!["CO".to_string()]);
        assert_eq!(preset.kind, ThermoScenarioKind::PlugFlowGas);
        assert!(preset.needs_thermochemistry);
        assert!(preset.needs_lambda);
        assert!(preset.needs_viscosity);
        assert!(preset.needs_diffusion);
        assert_eq!(
            preset.requested_properties(),
            vec![
                PropertyKind::Cp,
                PropertyKind::dH,
                PropertyKind::dS,
                PropertyKind::Lambda,
                PropertyKind::Visc,
            ]
        );
        assert_eq!(
            preset.output_policy,
            ThermoOutputPolicy::PropertiesAndTransport
        );
    }

    #[test]
    fn plug_flow_cond_preset_requests_conductivity_only() {
        let preset = ThermoScenarioPreset::plug_flow_cond(vec!["CO".to_string()]);
        assert_eq!(preset.kind, ThermoScenarioKind::PlugFlowCond);
        assert!(preset.needs_thermochemistry);
        assert!(preset.needs_lambda);
        assert!(!preset.needs_viscosity);
        assert!(!preset.needs_diffusion);
        assert_eq!(
            preset.output_policy,
            ThermoOutputPolicy::PropertiesAndTransport
        );
        assert_eq!(
            preset.requested_properties(),
            vec![
                PropertyKind::Cp,
                PropertyKind::dH,
                PropertyKind::dS,
                PropertyKind::Lambda,
            ]
        );
    }

    #[test]
    fn lookup_policy_and_preview_report_include_canonical_summary_rows() {
        let preset = ThermoScenarioPreset::plug_flow_cond(vec!["CO".to_string()])
            .with_priority_library(LibraryId::NasaGas)
            .with_permitted_library(LibraryId::Nist)
            .with_explicit_instruction("H2O".to_string(), LibraryId::NasaGas)
            .with_temperature_interval(300.0, 1500.0)
            .with_output_policy(ThermoOutputPolicy::PropertiesAndPreviewReport);

        let rows = preset.summary_rows();
        assert!(
            rows.iter()
                .any(|(key, value)| key == "scenario" && value == "PlugFlowCond")
        );
        assert!(
            rows.iter()
                .any(|(key, value)| key == "output_policy"
                    && value == "properties_and_preview_report")
        );
        assert!(
            rows.iter()
                .any(|(key, value)| key == "priority_libraries" && value == "NASA_gas")
        );
        assert!(
            rows.iter()
                .any(|(key, value)| key == "permitted_libraries" && value == "NIST")
        );
        assert!(
            rows.iter()
                .any(|(key, value)| key == "explicit_instructions"
                    && value.contains("H2O=>NASA_gas"))
        );
        assert!(
            rows.iter()
                .any(|(key, value)| key == "temperature_interval" && value == "300..1500")
        );

        let preview = preset.preview_report();
        assert_eq!(preview.kind, ThermoScenarioKind::PlugFlowCond);
        assert_eq!(
            preview.output_policy,
            ThermoOutputPolicy::PropertiesAndPreviewReport
        );
        assert!(
            preview
                .lines()
                .iter()
                .any(|line| line.contains("scenario: PlugFlowCond"))
        );
    }

    #[test]
    fn requested_inputs_are_stable_and_readable() {
        let preset = ThermoScenarioPreset::stirred_tank(vec!["CO".to_string()])
            .with_priority_library(LibraryId::NasaGas)
            .with_output_policy(ThermoOutputPolicy::PropertiesAndWriteWorkflow);

        let requested_inputs = preset.requested_inputs();
        assert!(
            requested_inputs
                .iter()
                .any(|line| line.contains("scenario: StirredTank"))
        );
        assert!(
            requested_inputs
                .iter()
                .any(|line| line.contains("output_policy: properties_and_write_workflow"))
        );
        assert!(
            requested_inputs
                .iter()
                .any(|line| line.contains("priority_libraries: NASA_gas"))
        );
    }

    #[test]
    fn multiple_substances_are_preserved_in_preview_and_data_request() {
        let preset = ThermoScenarioPreset::plug_flow_gas(vec![
            "CH4".to_string(),
            "O2".to_string(),
            "N2".to_string(),
        ])
        .with_priority_library(LibraryId::NasaGas)
        .with_permitted_library(LibraryId::Nist)
        .with_temperature(1200.0)
        .with_pressure(101325.0);

        let preview = preset.preview_report();
        let data_request = preset.requested_data_request();

        assert_eq!(
            preset.substances,
            vec!["CH4".to_string(), "O2".to_string(), "N2".to_string()]
        );
        assert!(
            preview
                .lines()
                .iter()
                .any(|line| line.contains("substances: CH4, O2, N2"))
        );
        assert_eq!(
            data_request.substances,
            vec!["CH4".to_string(), "O2".to_string(), "N2".to_string()]
        );
        assert_eq!(
            data_request.property_kinds,
            vec![
                PropertyKind::Cp,
                PropertyKind::dH,
                PropertyKind::dS,
                PropertyKind::Lambda,
                PropertyKind::Visc,
            ]
        );
        assert_eq!(
            data_request.numeric_data_types,
            vec![
                DataType::Cp,
                DataType::dH,
                DataType::dS,
                DataType::Lambda,
                DataType::Visc,
            ]
        );
        assert_eq!(
            data_request.function_data_types,
            vec![
                DataType::Cp_fun,
                DataType::dH_fun,
                DataType::dS_fun,
                DataType::Lambda_fun,
                DataType::Visc_fun,
            ]
        );
        assert_eq!(
            data_request.symbolic_data_types,
            vec![
                DataType::Cp_sym,
                DataType::dH_sym,
                DataType::dS_sym,
                DataType::Lambda_sym,
                DataType::Visc_sym,
            ]
        );
        assert!(data_request.requires_molar_masses);
    }

    #[test]
    fn lookup_policy_is_empty_by_default() {
        let policy = ThermoLookupPolicy::new();
        assert!(policy.is_empty());
        assert!(policy.priority_libraries.is_empty());
        assert!(policy.permitted_libraries.is_empty());
        assert!(policy.explicit_instructions.is_empty());
        assert!(policy.temperature_interval.is_none());
    }

    #[test]
    fn preset_can_apply_core_lookup_policy_without_forcing_thermo_range_lookup() {
        let mut subs_data = SubsData::new();
        let mut molar_masses = HashMap::new();
        molar_masses.insert("CO".to_string(), 28.0);

        let preset = ThermoScenarioPreset::stirred_tank(vec!["CO".to_string()])
            .with_priority_library(LibraryId::NasaGas)
            .with_permitted_library(LibraryId::Nist)
            .with_molar_masses(molar_masses)
            .with_temperature(298.15)
            .with_pressure(101325.0);

        preset.apply_to(&mut subs_data).unwrap();

        assert_eq!(subs_data.substances(), &["CO".to_string()]);
        assert_eq!(
            subs_data.library_priorities().get("NASA_gas"),
            Some(&LibraryPriority::Priority)
        );
        assert_eq!(
            subs_data.library_priorities().get("NIST"),
            Some(&LibraryPriority::Permitted)
        );
        assert_eq!(
            subs_data
                .explicit_search_instructions()
                .get("H2O")
                .map(String::as_str),
            None
        );
        assert_eq!(subs_data.temperature(), Some(298.15));
        assert_eq!(subs_data.pressure(), Some(101325.0));
        assert_eq!(subs_data.molar_mass_of("CO"), Some(28.0));
        assert_eq!(preset.lookup_policy.temperature_interval, None);
    }

    #[test]
    fn requested_data_request_includes_numeric_function_and_symbolic_layers() {
        let preset = ThermoScenarioPreset::plug_flow_cond(vec!["CO".to_string()]);
        let request = preset.requested_data_request();

        assert_eq!(
            request.numeric_data_types,
            vec![DataType::Cp, DataType::dH, DataType::dS, DataType::Lambda,]
        );
        assert_eq!(
            request.function_data_types,
            vec![
                DataType::Cp_fun,
                DataType::dH_fun,
                DataType::dS_fun,
                DataType::Lambda_fun,
            ]
        );
        assert_eq!(
            request.symbolic_data_types,
            vec![
                DataType::Cp_sym,
                DataType::dH_sym,
                DataType::dS_sym,
                DataType::Lambda_sym,
            ]
        );
    }
}
