//! Canonical physical-state declarations used by thermodynamic lookup.
//!
//! A physical state is a request constraint for a library record, not a
//! thermodynamic activity model.  For example, a liquid record may later be
//! evaluated as a pure condensed phase or as part of a solution model.

/// Requested or observed physical state of a thermodynamic record.
///
/// `Condensed` is an intentionally coarse request meaning either liquid or
/// solid.  It is useful for legacy callers but does not make an ambiguous
/// liquid/solid lookup deterministic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PhysicalState {
    Gas,
    Liquid,
    Solid,
    Condensed,
}

impl PhysicalState {
    /// Returns whether an observed record state satisfies this request.
    pub const fn accepts(self, observed: Self) -> bool {
        matches!(
            (self, observed),
            (Self::Gas, Self::Gas)
                | (Self::Liquid, Self::Liquid)
                | (Self::Solid, Self::Solid)
                | (
                    Self::Condensed,
                    Self::Liquid | Self::Solid | Self::Condensed
                )
                | (Self::Liquid | Self::Solid, Self::Condensed)
        )
    }
}

/// Typed lookup request.  The string remains the canonical substance name
/// supplied by an application; state is deliberately a separate field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermoRecordQuery {
    substance: String,
    physical_state: Option<PhysicalState>,
}

impl ThermoRecordQuery {
    pub fn new(substance: impl Into<String>) -> Self {
        Self {
            substance: substance.into(),
            physical_state: None,
        }
    }

    pub fn with_physical_state(mut self, physical_state: PhysicalState) -> Self {
        self.physical_state = Some(physical_state);
        self
    }

    pub fn with_physical_state_opt(mut self, physical_state: Option<PhysicalState>) -> Self {
        self.physical_state = physical_state;
        self
    }

    pub fn substance(&self) -> &str {
        &self.substance
    }

    pub const fn physical_state(&self) -> Option<PhysicalState> {
        self.physical_state
    }
}

/// Why a record was classified as belonging to a physical state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhysicalStateEvidence {
    /// The state is encoded in this library family's record key convention.
    KeyConvention,
    /// The library family represents one physical state unless a key says otherwise.
    LibraryDefault,
}
