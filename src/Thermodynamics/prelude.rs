//! Common thermodynamics imports for application-level code.
//!
//! This module is intentionally small and opinionated: it re-exports the
//! stable, frequently used pieces of the thermodynamics surface so callers can
//! write concise code without reaching into every submodule manually.

pub use crate::Thermodynamics::DBhandlers::thermo_api::{
    ThermoCalculator, ThermoEnum, ThermoError, create_thermal_by_name,
};
pub use crate::Thermodynamics::DBhandlers::transport_api::{
    TransportCalculator, TransportEnum, TransportError, create_transport_calculator_by_name,
};
pub use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, DerivedValueState, LibraryPriority, Phases, PropertyKind,
    PropertySearchState, SearchResult, SubsData, SubstanceSearchState, WhatIsFound,
};
pub use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
pub use crate::Thermodynamics::User_substances_presets::{
    ThermoLookupPolicy, ThermoOutputPolicy, ThermoPresetPreviewReport, ThermoRequestedData,
    ThermoScenarioKind, ThermoScenarioPreset,
};
pub use crate::Thermodynamics::physical_state::{
    PhysicalState, PhysicalStateEvidence, ThermoRecordQuery,
};
pub use crate::Thermodynamics::thermo_lib_api::{
    LibraryCapability, LibraryId, ResolvedThermoRecord, ThermoData, ThermoLibraryError,
    ThermoRepository,
};

#[cfg(test)]
mod tests {
    #[test]
    fn prelude_exports_core_thermodynamics_types() {
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::SubsData>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::ThermoData>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::LibraryId>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::ThermoLookupPolicy>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::ThermoOutputPolicy>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::ThermoPresetPreviewReport>();
        let _ = std::any::type_name::<crate::Thermodynamics::prelude::ThermoRequestedData>();
    }
}
