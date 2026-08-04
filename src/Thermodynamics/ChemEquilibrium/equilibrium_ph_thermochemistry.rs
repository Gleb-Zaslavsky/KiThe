//! Typed thermochemistry capabilities for the fixed-pressure, fixed-enthalpy engine.
//!
//! This module is intentionally below workflow orchestration. It resolves
//! component-aligned Gibbs, enthalpy, heat-capacity, interval, and provenance
//! capabilities from already selected phase data, but it does not own solver
//! policy, phase lifecycle, scalar bracketing, or result publication.

use std::fmt;
use std::rc::Rc;
use std::sync::{Arc, Mutex};

use crate::Thermodynamics::phase_layout::PhaseComponentId;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::User_PhaseOrSolution::ResolvedPhaseSystem;
use crate::Thermodynamics::User_substances::{CalculatorType, DataType, SubsData, WhatIsFound};
use RustedSciThe::symbolic::symbolic_engine::Expr;

/// One temperature-dependent molar enthalpy capability in J/mol.
pub type MolarEnthalpyFunction<'a> =
    Arc<dyn Fn(f64) -> Result<f64, ReactionExtentError> + Send + Sync + 'a>;

/// One temperature-dependent molar thermochemistry value in J/mol (or J/mol/K
/// for heat capacity). The fallible return keeps record-domain failures
/// visible to the outer `P,H` workflow instead of turning them into NaNs.
pub type MolarThermoFunction =
    Arc<dyn Fn(f64) -> Result<f64, ReactionExtentError> + Send + Sync + 'static>;

/// Provenance for one phase-qualified thermochemical component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermochemistryProvenance {
    component: PhaseComponentId,
    library: String,
    record_key: String,
    state: String,
}

/// Symbolic standard-state capabilities in the exact resolved component order.
///
/// Numeric closures remain the canonical, format-agnostic thermochemistry
/// interface. This optional companion exists only for backends, such as RST,
/// that can construct and differentiate a complete nonlinear system from
/// expressions. It is prepared from the same private `SubsData` copy as the
/// numeric capabilities, so it cannot silently select different records.
///
/// A bundle represents one exact native polynomial interval for every
/// component. NASA and NIST records may change coefficients at a temperature
/// boundary; the symbolic P,H path must never extend one interval's expression
/// across that boundary and pretend that it is globally valid.
#[derive(Clone, Debug)]
pub(crate) struct ResolvedSymbolicThermochemistry {
    standard_gibbs: Vec<Expr>,
    enthalpy: Vec<Expr>,
}

impl ResolvedSymbolicThermochemistry {
    pub(crate) fn new(
        standard_gibbs: Vec<Expr>,
        enthalpy: Vec<Expr>,
    ) -> Result<Self, ReactionExtentError> {
        if standard_gibbs.is_empty() || standard_gibbs.len() != enthalpy.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic thermochemistry has {} Gibbs and {} enthalpy expressions",
                standard_gibbs.len(),
                enthalpy.len()
            )));
        }
        Ok(Self {
            standard_gibbs,
            enthalpy,
        })
    }

    pub(crate) fn standard_gibbs(&self) -> &[Expr] {
        &self.standard_gibbs
    }

    pub(crate) fn enthalpy(&self) -> &[Expr] {
        &self.enthalpy
    }

    /// Extracts a reduced symbolic thermochemistry bundle for an active subset.
    ///
    /// Returns a new [`PhSymbolicThermochemistry`] containing only the Gibbs
    /// and enthalpy expressions at the requested indices. This is used by the
    /// active-set projection path where the monolithic `P,H` formulation must
    /// rebuild its symbolic closures for a reduced phase set without touching
    /// the repository or re-parsing thermochemical records.
    fn subset(&self, indices: &[usize]) -> Result<Self, ReactionExtentError> {
        let mut standard_gibbs = Vec::with_capacity(indices.len());
        let mut enthalpy = Vec::with_capacity(indices.len());
        for &index in indices {
            let gibbs = self.standard_gibbs.get(index).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "symbolic P,H active-set index {index} exceeds Gibbs capability count {}",
                    self.standard_gibbs.len()
                ))
            })?;
            let heat = self.enthalpy.get(index).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "symbolic P,H active-set index {index} exceeds enthalpy capability count {}",
                    self.enthalpy.len()
                ))
            })?;
            standard_gibbs.push(gibbs.clone());
            enthalpy.push(heat.clone());
        }
        Self::new(standard_gibbs, enthalpy)
    }
}

impl ThermochemistryProvenance {
    /// Creates provenance for a component-aligned capability supplied by an
    /// adapter or a synthetic test fixture.
    pub fn new(
        component: PhaseComponentId,
        library: impl Into<String>,
        record_key: impl Into<String>,
        state: impl Into<String>,
    ) -> Self {
        Self {
            component,
            library: library.into(),
            record_key: record_key.into(),
            state: state.into(),
        }
    }

    /// Phase-qualified component identity.
    pub fn component(&self) -> &PhaseComponentId {
        &self.component
    }

    /// Selected thermochemistry library.
    pub fn library(&self) -> &str {
        &self.library
    }

    /// Exact record key selected in that library.
    pub fn record_key(&self) -> &str {
        &self.record_key
    }

    /// State evidence captured by the lookup report.
    pub fn state(&self) -> &str {
        &self.state
    }
}

/// Thermochemistry capabilities and provenance for one solver layout.
///
/// All vectors use exactly the same order as `ResolvedPhaseSystem::layout()`.
/// `g0(T)` and `h(T)` are mandatory; `Cp(T)` is optional until an analytic
/// temperature derivative or a monolithic `P,H` formulation needs it. The
/// bundle is deliberately independent of the source polynomial format: its
/// input is the capability API of the already resolved `ThermoCalculator`,
/// not a NASA/NIST coefficient representation.
#[derive(Clone)]
pub struct ResolvedThermochemistry {
    provenance: Vec<ThermochemistryProvenance>,
    temperature_bounds: TemperatureBounds,
    gibbs: Vec<MolarThermoFunction>,
    enthalpy: Vec<MolarThermoFunction>,
    heat_capacity: Vec<Option<MolarThermoFunction>>,
    /// Explicit symbolic expressions installed by an adapter or synthetic
    /// fixture. Real `SubsData` expressions are built lazily for the actual
    /// P,H bounds from `symbolic_sources` below.
    symbolic: Option<ResolvedSymbolicThermochemistry>,
    /// Private phase-local data copies used only to materialise an exact
    /// single-interval symbolic payload for an RST P,H solve. They are kept
    /// separate from numeric closure caches so symbolic preparation cannot
    /// mutate canonical numeric evaluation state.
    symbolic_sources: Option<std::collections::HashMap<Option<String>, SubsData>>,
}

impl fmt::Debug for ResolvedThermochemistry {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ResolvedThermochemistry")
            .field("provenance", &self.provenance)
            .field("temperature_bounds", &self.temperature_bounds)
            .field("component_count", &self.len())
            .finish_non_exhaustive()
    }
}

impl PartialEq for ResolvedThermochemistry {
    fn eq(&self, other: &Self) -> bool {
        self.provenance == other.provenance && self.temperature_bounds == other.temperature_bounds
    }
}

impl ResolvedThermochemistry {
    /// Builds a bundle from the records already selected by phase resolution.
    ///
    /// Each phase is evaluated in a private SubsData clone. Gibbs, enthalpy,
    /// heat-capacity, interval, and provenance therefore all refer to the same
    /// selected record without mutating the resolved system or repository.
    /// The method consumes the unified calculator capability API and does not
    /// inspect a NASA/NIST polynomial layout.
    pub fn from_resolved_system(
        resolved: &ResolvedPhaseSystem,
    ) -> Result<Self, ReactionExtentError> {
        struct PhaseSource {
            data: Arc<Mutex<SubsData>>,
            gibbs: Arc<Mutex<PhaseGibbsSource>>,
            provenance: std::collections::HashMap<String, ThermochemistryProvenance>,
        }

        let mut phase_sources = std::collections::HashMap::new();
        let mut symbolic_sources = std::collections::HashMap::new();
        let mut common_bounds: Option<TemperatureBounds> = None;

        for (phase, original) in resolved.phase_data() {
            let local = original.clone();
            let data = Arc::new(Mutex::new(local));
            let gibbs = Arc::new(Mutex::new(PhaseGibbsSource::new(original.clone())));
            let substances = original.substances().to_vec();
            let mut phase_bounds: Option<TemperatureBounds> = None;
            let mut provenance = std::collections::HashMap::new();

            for substance in substances {
                let source = data
                    .lock()
                    .map_err(|_| ReactionExtentError::InvalidProblem {
                        field: "thermochemistry_source",
                        message: "private SubsData source was poisoned during preparation"
                            .to_string(),
                    })?;
                let record = source
                    .get_search_result(&substance, WhatIsFound::Thermo)
                    .ok_or_else(|| ReactionExtentError::InvalidProblem {
                        field: "thermochemistry_provenance",
                        message: format!("no thermochemistry record for '{substance}'"),
                    })?;
                let mut calculator = match record.calculator().cloned() {
                    Some(CalculatorType::Thermo(calculator)) => calculator,
                    Some(CalculatorType::Transport(_)) => {
                        return Err(ReactionExtentError::InvalidProblem {
                            field: "thermochemistry_interval",
                            message: format!(
                                "transport calculator cannot provide thermochemistry for '{substance}'"
                            ),
                        });
                    }
                    None => {
                        return Err(ReactionExtentError::InvalidProblem {
                            field: "thermochemistry_interval",
                            message: format!("calculator is missing for '{substance}'"),
                        });
                    }
                };
                calculator.parse_coefficients().map_err(|error| {
                    ReactionExtentError::InvalidProblem {
                        field: "thermochemistry_interval",
                        message: format!(
                            "failed to parse temperature domain for '{substance}': {error}"
                        ),
                    }
                })?;
                let interval = calculator.valid_temperature_interval().map_err(|error| {
                    ReactionExtentError::InvalidProblem {
                        field: "thermochemistry_interval",
                        message: format!(
                            "failed to read temperature domain for '{substance}': {error}"
                        ),
                    }
                })?;
                let interval = TemperatureBounds::new(interval.0, interval.1)?;
                phase_bounds = Some(match phase_bounds {
                    Some(current) => current.intersect(interval).map_err(|_| {
                        ReactionExtentError::InvalidProblem {
                            field: "thermochemistry_interval",
                            message: format!(
                                "temperature domains do not overlap in phase {phase:?}"
                            ),
                        }
                    })?,
                    None => interval,
                });

                let state = resolved
                    .report()
                    .phase(&crate::Thermodynamics::phase_layout::PhaseId::new(
                        phase.clone(),
                    ))
                    .and_then(|summary| {
                        summary.search().rows().iter().find(|row| {
                            row.substance() == substance
                                && row.library() == record.library()
                                && row.record_key() == record.record_key()
                        })
                    })
                    .map(|row| row.state().to_string())
                    .unwrap_or_else(|| "resolved".to_string());
                provenance.insert(
                    substance.clone(),
                    ThermochemistryProvenance::new(
                        PhaseComponentId::new(
                            crate::Thermodynamics::phase_layout::PhaseId::new(phase.clone()),
                            substance.clone(),
                        ),
                        record.library(),
                        record.record_key(),
                        state,
                    ),
                );
            }

            let bounds = phase_bounds.ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_interval",
                message: format!("phase {phase:?} has no thermochemistry components"),
            })?;
            common_bounds =
                Some(match common_bounds {
                    Some(current) => current.intersect(bounds).map_err(|_| {
                        ReactionExtentError::InvalidProblem {
                            field: "thermochemistry_interval",
                            message: "selected records have no common temperature domain"
                                .to_string(),
                        }
                    })?,
                    None => bounds,
                });
            // Symbolic expressions are materialised only after the P,H bounds
            // are known. A NASA/NIST record can change its native polynomial
            // at an internal temperature boundary, so building one expression
            // here would give it an unjustified global validity range.
            symbolic_sources.insert(phase.clone(), original.clone());
            phase_sources.insert(
                phase.clone(),
                PhaseSource {
                    data,
                    gibbs,
                    provenance,
                },
            );
        }

        let count = resolved.layout().component_count();
        let mut provenance = Vec::with_capacity(count);
        let mut gibbs = Vec::with_capacity(count);
        let mut enthalpy = Vec::with_capacity(count);
        let mut heat_capacity = Vec::with_capacity(count);
        for component in resolved.layout().components() {
            let phase = phase_sources
                .get_mut(component.phase.as_option())
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "thermochemistry_layout",
                    message: format!("missing thermochemistry phase for '{}'", component.label()),
                })?;
            let row = phase
                .provenance
                .get(&component.substance)
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "thermochemistry_provenance",
                    message: format!("missing provenance for '{}'", component.label()),
                })?
                .clone();
            provenance.push(row);
            let source = Arc::clone(&phase.data);
            let substance = component.substance.clone();
            gibbs.push(gibbs_function(Arc::clone(&phase.gibbs), substance.clone()));
            enthalpy.push(property_function(
                Arc::clone(&source),
                substance.clone(),
                DataType::dH_fun,
            ));
            heat_capacity.push(Some(property_function(source, substance, DataType::Cp_fun)));
        }

        let mut thermochemistry = Self::from_functions(
            provenance,
            common_bounds.ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_interval",
                message: "resolved system has no thermochemistry domain".to_string(),
            })?,
            gibbs,
            enthalpy,
            heat_capacity,
        )?;
        thermochemistry.symbolic_sources = Some(symbolic_sources);
        Ok(thermochemistry)
    }

    /// Builds a validated bundle from component-aligned capabilities.
    pub fn from_functions(
        provenance: Vec<ThermochemistryProvenance>,
        temperature_bounds: TemperatureBounds,
        gibbs: Vec<MolarThermoFunction>,
        enthalpy: Vec<MolarThermoFunction>,
        heat_capacity: Vec<Option<MolarThermoFunction>>,
    ) -> Result<Self, ReactionExtentError> {
        let count = provenance.len();
        if count == 0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "thermochemistry_bundle",
                message: "at least one thermochemical component is required".to_string(),
            });
        }
        if gibbs.len() != count || enthalpy.len() != count || heat_capacity.len() != count {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "thermochemistry bundle has {} provenance rows, {} Gibbs, {} enthalpy, and {} Cp capabilities",
                count,
                gibbs.len(),
                enthalpy.len(),
                heat_capacity.len()
            )));
        }
        Ok(Self {
            provenance,
            temperature_bounds,
            gibbs,
            enthalpy,
            heat_capacity,
            symbolic: None,
            symbolic_sources: None,
        })
    }

    /// Number of phase-qualified components.
    pub fn len(&self) -> usize {
        self.provenance.len()
    }

    /// Whether this bundle contains no components.
    pub fn is_empty(&self) -> bool {
        self.provenance.is_empty()
    }

    /// Common valid temperature interval for all selected records.
    pub fn temperature_bounds(&self) -> TemperatureBounds {
        self.temperature_bounds
    }

    /// Component provenance in solver order.
    pub fn provenance(&self) -> &[ThermochemistryProvenance] {
        &self.provenance
    }

    /// Creates a component-aligned view for one active-set projection.
    ///
    /// The projection owns the canonical active species order, so preserving
    /// that order here keeps Gibbs, enthalpy, Cp, and provenance rows aligned
    /// with the reduced nonlinear problem.
    pub(crate) fn subset(&self, indices: &[usize]) -> Result<Self, ReactionExtentError> {
        let mut provenance = Vec::with_capacity(indices.len());
        let mut gibbs = Vec::with_capacity(indices.len());
        let mut enthalpy = Vec::with_capacity(indices.len());
        let mut heat_capacity = Vec::with_capacity(indices.len());
        for &index in indices {
            let row = self.provenance.get(index).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "P,H active-set index {index} exceeds thermochemistry component count {}",
                    self.len()
                ))
            })?;
            provenance.push(row.clone());
            gibbs.push(self.gibbs[index].clone());
            enthalpy.push(self.enthalpy[index].clone());
            heat_capacity.push(self.heat_capacity[index].clone());
        }
        let mut subset = Self::from_functions(
            provenance,
            self.temperature_bounds,
            gibbs,
            enthalpy,
            heat_capacity,
        )?;
        subset.symbolic = self
            .symbolic
            .as_ref()
            .map(|symbolic| symbolic.subset(indices))
            .transpose()?;
        subset.symbolic_sources = self.symbolic_sources.clone();
        Ok(subset)
    }

    /// Builds exact symbolic thermochemistry for one bounded P,H solve.
    ///
    /// Real NASA/NIST sources are pinned to the lower bound first. The upper
    /// bound must remain in the same native coefficient interval for every
    /// component. This keeps RST's residual/Jacobian equivalent to the
    /// canonical piecewise numeric route. Crossing a boundary is reported as
    /// an unavailable symbolic capability rather than fitting or extrapolating
    /// behind the caller's back.
    pub(crate) fn symbolic_for_bounds(
        &self,
        bounds: TemperatureBounds,
    ) -> Result<Option<ResolvedSymbolicThermochemistry>, String> {
        if let Some(symbolic) = &self.symbolic {
            return Ok(Some(symbolic.clone()));
        }
        let Some(sources) = &self.symbolic_sources else {
            return Ok(None);
        };

        let mut substances_by_phase: std::collections::HashMap<Option<String>, Vec<String>> =
            std::collections::HashMap::new();
        for row in &self.provenance {
            substances_by_phase
                .entry(row.component().phase.as_option().clone())
                .or_default()
                .push(row.component().substance.clone());
        }

        let mut expressions_by_phase = std::collections::HashMap::new();
        for (phase, source) in sources {
            let substances = substances_by_phase
                .get(phase)
                .ok_or_else(|| format!("missing symbolic component selection for '{phase:?}'"))?;
            let expressions = symbolic_phase_capabilities(source, substances, bounds)?;
            expressions_by_phase.insert(phase.clone(), expressions);
        }

        let mut standard_gibbs = Vec::with_capacity(self.len());
        let mut enthalpy = Vec::with_capacity(self.len());
        for row in &self.provenance {
            let phase = row.component().phase.as_option();
            let substance = &row.component().substance;
            let (gibbs_by_substance, enthalpy_by_substance) = expressions_by_phase
                .get(phase)
                .ok_or_else(|| format!("missing symbolic phase source for '{phase:?}'"))?;
            let gibbs = gibbs_by_substance.get(substance).ok_or_else(|| {
                format!(
                    "missing symbolic standard Gibbs expression for '{}'",
                    row.component().label()
                )
            })?;
            let heat = enthalpy_by_substance.get(substance).ok_or_else(|| {
                format!(
                    "missing symbolic enthalpy expression for '{}'",
                    row.component().label()
                )
            })?;
            standard_gibbs.push(gibbs.clone());
            enthalpy.push(heat.clone());
        }
        ResolvedSymbolicThermochemistry::new(standard_gibbs, enthalpy)
            .map(Some)
            .map_err(|error| error.to_string())
    }

    #[cfg(test)]
    pub(crate) fn with_symbolic_expressions(
        mut self,
        standard_gibbs: Vec<Expr>,
        enthalpy: Vec<Expr>,
    ) -> Result<Self, ReactionExtentError> {
        if standard_gibbs.len() != self.len() || enthalpy.len() != self.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "symbolic test capability has {} Gibbs and {} enthalpy expressions for {} components",
                standard_gibbs.len(),
                enthalpy.len(),
                self.len()
            )));
        }
        self.symbolic = Some(ResolvedSymbolicThermochemistry::new(
            standard_gibbs,
            enthalpy,
        )?);
        Ok(self)
    }

    /// Evaluates standard Gibbs energies in solver order.
    pub fn evaluate_gibbs(&self, temperature: f64) -> Result<Vec<f64>, ReactionExtentError> {
        self.evaluate_functions("standard_gibbs", &self.gibbs, temperature)
    }

    /// Materializes a finite Gibbs snapshot for the retained infallible
    /// phase-stability boundary.
    ///
    /// The old phase-control helper accepts `GibbsFn = Fn(f64) -> f64`, while
    /// the resolved P,H bundle deliberately exposes
    /// `Fn(f64) -> Result<f64, ReactionExtentError>`.  Converting a fallible
    /// capability with `unwrap_or(NAN)` would hide a broken thermochemical
    /// record until much later in the nonlinear solve.  Instead, evaluate all
    /// values first and cross the compatibility boundary only after this
    /// method has validated the complete snapshot.
    ///
    /// The returned closures are intentionally constant.  They are suitable
    /// only for consumers that evaluate phase stability at the same accepted
    /// temperature supplied here; callers must not use them as a general
    /// temperature-dependent thermochemistry API.
    pub(crate) fn gibbs_snapshot_for_legacy_boundary(
        &self,
        temperature: f64,
    ) -> Result<Vec<GibbsFn>, ReactionExtentError> {
        let values = self.evaluate_gibbs(temperature)?;
        Ok(values
            .into_iter()
            .map(|value| Rc::new(move |_| value) as GibbsFn)
            .collect())
    }

    /// Evaluates molar enthalpies in solver order.
    pub fn evaluate_enthalpy(&self, temperature: f64) -> Result<Vec<f64>, ReactionExtentError> {
        self.evaluate_functions("molar_enthalpy", &self.enthalpy, temperature)
    }

    /// Evaluates optional heat capacities in solver order.
    pub fn evaluate_heat_capacity(
        &self,
        temperature: f64,
    ) -> Result<Vec<Option<f64>>, ReactionExtentError> {
        self.temperature_bounds
            .contains(temperature)
            .then_some(())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "temperature",
                message: format!(
                    "temperature {temperature} K is outside the common thermochemistry interval"
                ),
            })?;
        self.heat_capacity
            .iter()
            .enumerate()
            .map(|(index, function)| {
                function
                    .as_ref()
                    .map(|function| {
                        let value = function(temperature)?;
                        if !value.is_finite() {
                            return Err(ReactionExtentError::InvalidProblem {
                                field: "heat_capacity",
                                message: format!("Cp function {index} returned a non-finite value"),
                            });
                        }
                        Ok(value)
                    })
                    .transpose()
            })
            .collect()
    }

    /// Returns the mandatory enthalpy capabilities for the compatibility
    /// enthalpy model owned by the workflow facade.
    pub(crate) fn enthalpy_functions(&self) -> Vec<MolarEnthalpyFunction<'static>> {
        self.enthalpy.clone()
    }

    /// Returns the optional Cp capabilities for the compatibility enthalpy
    /// model owned by the workflow facade.
    pub(crate) fn heat_capacity_functions(&self) -> Vec<Option<MolarEnthalpyFunction<'static>>> {
        self.heat_capacity.clone()
    }

    /// Evaluates one family of temperature-dependent molar functions for all
    /// components, returning values in solver order.
    ///
    /// The `field` name is used only for diagnostic error messages when a
    /// function returns a non-finite or out-of-bounds result. The temperature
    /// is validated against the common thermochemistry interval before any
    /// function is called, so callers receive a single typed error rather than
    /// a per-component panic or silent `NaN` propagation.
    fn evaluate_functions(
        &self,
        field: &'static str,
        functions: &[MolarThermoFunction],
        temperature: f64,
    ) -> Result<Vec<f64>, ReactionExtentError> {
        self.temperature_bounds
            .contains(temperature)
            .then_some(())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "temperature",
                message: format!(
                    "temperature {temperature} K is outside the common thermochemistry interval"
                ),
            })?;
        functions
            .iter()
            .enumerate()
            .map(|(index, function)| {
                let value = function(temperature)?;
                if !value.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field,
                        message: format!(
                            "thermochemistry function {index} returned a non-finite value"
                        ),
                    });
                }
                Ok(value)
            })
            .collect()
    }
}

fn property_function(
    data: Arc<Mutex<SubsData>>,
    substance: String,
    data_type: DataType,
) -> MolarThermoFunction {
    Arc::new(move |temperature| {
        let mut data = data
            .lock()
            .map_err(|_| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_source",
                message: "private SubsData source was poisoned during evaluation".to_string(),
            })?;
        // Do not use the copied closure cache here. A cloned SubsData may
        // intentionally contain shape-preserving zero placeholders, while
        // coefficient selection is temperature-dependent. The calculator is
        // the authoritative read path for this private source snapshot.
        data.extract_thermal_coeffs(&substance, temperature)
            .map_err(|error| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_property",
                message: format!("failed to select coefficients for '{substance}': {error}"),
            })?;
        let (cp, dh, _) = data
            .calculate_thermo_properties(&substance, temperature)
            .map_err(|error| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_property",
                message: format!("failed to evaluate '{substance}': {error}"),
            })?;
        match data_type {
            DataType::Cp_fun => Ok(cp),
            DataType::dH_fun => Ok(dh),
            _ => Err(ReactionExtentError::InvalidProblem {
                field: "thermochemistry_property",
                message: "unsupported property requested by the P,H bundle".to_string(),
            }),
        }
    })
}

/// Materializes symbolic Gibbs and enthalpy expressions from one private phase
/// copy. The caller owns this copy, so symbolic preparation cannot publish a
/// cache into the resolved system or alter a JSON-backed repository.
fn symbolic_phase_capabilities(
    original: &SubsData,
    selected_substances: &[String],
    bounds: TemperatureBounds,
) -> Result<
    (
        std::collections::HashMap<String, Expr>,
        std::collections::HashMap<String, Expr>,
    ),
    String,
> {
    let mut working = original.clone();
    // Preserve the resolved search/calculator maps in the clone, but restrict
    // all batch APIs to the active projection. An inactive component must not
    // deny a valid symbolic P,H solve merely because its own polynomial range
    // differs from the active system's range.
    working.substances = selected_substances.to_vec();
    working
        .extract_all_thermal_coeffs(bounds.lower())
        .map_err(|error| error.to_string())?;
    for substance in working.substances().to_vec() {
        let remains_in_native_interval = working
            .is_coeffs_valid_for_T(&substance, bounds.upper())
            .map_err(|error| error.to_string())?;
        if !remains_in_native_interval {
            return Err(format!(
                "symbolic P,H interval [{}, {}] crosses a native thermochemical coefficient boundary for '{substance}'",
                bounds.lower(),
                bounds.upper()
            ));
        }
    }
    working
        .calculate_therm_map_of_sym()
        .map_err(|error| error.to_string())?;
    let standard_gibbs = working
        .calculate_dG0_sym_one_phase()
        .map_err(|error| error.to_string())?;
    let mut enthalpy = std::collections::HashMap::new();
    for substance in working.substances().iter() {
        let expression = working
            .get_thermo_symbolic(substance, DataType::dH_sym)
            .ok_or_else(|| format!("missing symbolic dH expression for '{substance}'"))?;
        enthalpy.insert(substance.clone(), expression.clone());
    }
    Ok((standard_gibbs, enthalpy))
}

/// Per-phase standard-Gibbs source used only by the resolved P,H bridge.
///
/// A `SubsData` Gibbs closure snapshots selected polynomial coefficients. This
/// cache keeps one such snapshot for all components of a phase at a given
/// temperature, so a residual/Jacobian evaluation does not rebuild identical
/// closures once per component. The source is private to the bundle and never
/// publishes coefficient changes back into the resolved phase system. A new
/// temperature first refreshes all selected intervals, keeping P,H Gibbs data
/// identical to the canonical P,T bridge at polynomial boundaries.
struct PhaseGibbsSource {
    data: SubsData,
    selected_temperature: Option<f64>,
    functions: std::collections::HashMap<String, Box<dyn Fn(f64) -> f64 + Send + Sync>>,
}

impl PhaseGibbsSource {
    fn new(data: SubsData) -> Self {
        Self {
            data,
            selected_temperature: None,
            functions: std::collections::HashMap::new(),
        }
    }

    fn evaluate(&mut self, substance: &str, temperature: f64) -> Result<f64, ReactionExtentError> {
        if self.selected_temperature != Some(temperature) {
            self.data
                .extract_all_thermal_coeffs(temperature)
                .map_err(|error| ReactionExtentError::InvalidProblem {
                    field: "thermochemistry_gibbs",
                    message: format!(
                        "failed to select phase coefficients at {temperature} K: {error}"
                    ),
                })?;
            self.functions = self.data.calculate_dG0_fun_one_phase().map_err(|error| {
                ReactionExtentError::InvalidProblem {
                    field: "thermochemistry_gibbs",
                    message: format!(
                        "failed to build phase standard Gibbs functions at {temperature} K: {error}"
                    ),
                }
            })?;
            self.selected_temperature = Some(temperature);
        }
        let function =
            self.functions
                .get(substance)
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "thermochemistry_gibbs",
                    message: format!("missing standard Gibbs function for '{substance}'"),
                })?;
        let value = function(temperature);
        if !value.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "thermochemistry_gibbs",
                message: format!(
                    "standard Gibbs function for '{substance}' returned a non-finite value"
                ),
            });
        }
        Ok(value)
    }
}

fn gibbs_function(data: Arc<Mutex<PhaseGibbsSource>>, substance: String) -> MolarThermoFunction {
    Arc::new(move |temperature| {
        let mut data = data
            .lock()
            .map_err(|_| ReactionExtentError::InvalidProblem {
                field: "thermochemistry_source",
                message: "private SubsData source was poisoned during evaluation".to_string(),
            })?;
        data.evaluate(&substance, temperature)
    })
}
