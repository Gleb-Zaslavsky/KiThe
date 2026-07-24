//! Symbolic layout and revision-aligned cache views for phase thermodynamics.
//!
//! These types describe derived state only.  `PhaseSystem` owns mutation and
//! invalidation; facades expose the borrowed snapshots defined here.
//!
//! ```text
//! IERARCHY:
//!
//!                PhaseThermoCacheBundle  - Cache bundle for one phase. Gibbs and entropy are independent property families, each with aligned numeric, function, and symbolic views.
//!                /             \
//!           GIBBS             ENTROPY
//!             |                  |
//!        PhasePropertyCache      PhasePropertyCache  
//!        /      |       \          /          |    \   
//!   NUMERIC  FUNCTION  SYMBOLIC   NUMERIC  FUNCTION  SYMBOLIC
//!               |                             |
//!        PhaseFunction                   PhaseFunction    
//!
//! PhaseSymbolicLayout - Canonical symbolic mole-variable snapshot shared by every phase facade.
//!        PhaseSymbolicLayout
//!             |
//!        SystemLayout
//!            /      |         \
//!          phases   components  phase_ranges
//!          /         |             \
//!                 PhaseComponentId
//! ```

use std::collections::HashMap;
use std::collections::hash_map;
use std::sync::Arc;

use RustedSciThe::symbolic::symbolic_engine::Expr;

use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use crate::Thermodynamics::phase_layout::SystemLayout;

use super::{IndexedMoleVariables, PhaseFunction};

pub(crate) fn boxify_fun_map(
    fun_map: &HashMap<String, PhaseFunction>,
) -> HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>> {
    fun_map
        .iter()
        .map(|(name, fun)| {
            let fun = Arc::clone(fun);
            (
                name.clone(),
                Box::new(move |T: f64, x: Option<Vec<f64>>, P: Option<f64>| fun(T, x, P))
                    as Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>,
            )
        })
        .collect()
}

/// Canonical symbolic mole-variable snapshot shared by every phase facade.
///
/// Variable names are derived data, while `PhaseSystem` remains the sole owner
/// of the canonical component order that generated them. Keeping this snapshot
/// layout-free prevents two independently mutable sources of ordering truth.
#[derive(Clone, Debug, Default)]
pub(crate) struct PhaseSymbolicLayout {
    /// Phase-level variables: (phase mole fraction, [component mole fractions])
    phase_variables: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    /// Component-level variables: [component mole fractions]
    component_variables: Vec<Expr>,
    /// Phase totals: [phase mole fractions]
    phase_totals: Vec<Expr>,
    /// Variables indexed by phase and substance: {phase_name: {substance_name: variable}}
    variables_by_phase_and_substance: HashMap<Option<String>, HashMap<String, Expr>>,
}

impl PhaseSymbolicLayout {
    pub(crate) fn new(
        phase_variables: HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
        component_variables: Vec<Expr>,
        phase_totals: Vec<Expr>,
        variables_by_phase_and_substance: HashMap<Option<String>, HashMap<String, Expr>>,
    ) -> Self {
        Self {
            phase_variables,
            component_variables,
            phase_totals,
            variables_by_phase_and_substance,
        }
    }

    pub(crate) fn phase_variables(
        &self,
    ) -> &HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)> {
        &self.phase_variables
    }

    #[cfg(test)]
    pub(crate) fn component_variables(&self) -> &[Expr] {
        &self.component_variables
    }

    #[cfg(test)]
    pub(crate) fn phase_totals(&self) -> &[Expr] {
        &self.phase_totals
    }

    #[cfg(test)]
    pub(crate) fn variables_by_phase_and_substance(
        &self,
    ) -> &HashMap<Option<String>, HashMap<String, Expr>> {
        &self.variables_by_phase_and_substance
    }

    #[cfg(test)]
    pub(crate) fn is_empty(&self) -> bool {
        self.phase_variables.is_empty()
            && self.component_variables.is_empty()
            && self.phase_totals.is_empty()
            && self.variables_by_phase_and_substance.is_empty()
    }

    pub(crate) fn legacy_projection(&self) -> IndexedMoleVariables {
        (
            self.phase_variables.clone(),
            self.component_variables.clone(),
            self.phase_totals.clone(),
            self.variables_by_phase_and_substance.clone(),
        )
    }
}

/// All cached representations of one thermodynamic property in one phase.
///
/// Numeric values, executable functions, and symbolic expressions are three
/// views of the same property family. Keeping them together prevents the
/// engine from treating Gibbs and entropy caches as six unrelated maps.
#[derive(Clone, Default)]
pub(crate) struct PhasePropertyCache {
    pub(crate) numeric: HashMap<String, f64>,
    pub(crate) functions: HashMap<String, PhaseFunction>,
    pub(crate) symbolic: HashMap<String, Expr>,
}

impl std::fmt::Debug for PhasePropertyCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PhasePropertyCache")
            .field("numeric_len", &self.numeric.len())
            .field("function_len", &self.functions.len())
            .field("symbolic_len", &self.symbolic.len())
            .finish()
    }
}

/// Cache bundle for one phase. Gibbs and entropy are independent property
/// families, each with aligned numeric, function, and symbolic views.
#[derive(Clone, Default)]
pub(crate) struct PhaseThermoCacheBundle {
    pub(crate) gibbs: PhasePropertyCache,
    pub(crate) entropy: PhasePropertyCache,
}

impl PhaseThermoCacheBundle {
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

impl std::fmt::Debug for PhaseThermoCacheBundle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PhaseThermoCacheBundle")
            .field("gibbs", &self.gibbs)
            .field("entropy", &self.entropy)
            .finish()
    }
}

/// Borrowed view over a property projected from nested phase caches.
#[derive(Clone, Copy)]
pub struct PhaseThermoPropertyView<'a, T> {
    phase_bundles: &'a HashMap<Option<String>, PhaseThermoCacheBundle>,
    selector: fn(&PhaseThermoCacheBundle) -> &HashMap<String, T>,
}

/// Allocation-free iterator over one property cache for every canonical phase.
/// It borrows phase identities and property maps directly from their bundles.
pub struct PhaseThermoPropertyIter<'a, T> {
    inner: hash_map::Iter<'a, Option<String>, PhaseThermoCacheBundle>,
    selector: fn(&PhaseThermoCacheBundle) -> &HashMap<String, T>,
}

impl<'a, T: 'a> Iterator for PhaseThermoPropertyIter<'a, T> {
    type Item = (&'a Option<String>, &'a HashMap<String, T>);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|(phase, bundle)| (phase, (self.selector)(bundle)))
    }
}

impl<'a, T> PhaseThermoPropertyView<'a, T> {
    pub(crate) fn new(
        phase_bundles: &'a HashMap<Option<String>, PhaseThermoCacheBundle>,
        selector: fn(&PhaseThermoCacheBundle) -> &HashMap<String, T>,
    ) -> Self {
        Self {
            phase_bundles,
            selector,
        }
    }

    pub fn get(&self, phase: &Option<String>) -> Option<&HashMap<String, T>> {
        self.phase_bundles
            .get(phase)
            .map(|bundle| (self.selector)(bundle))
    }

    pub fn contains_key(&self, phase: &Option<String>) -> bool {
        self.get(phase).is_some()
    }

    pub fn is_empty(&self) -> bool {
        self.phase_bundles
            .values()
            .all(|bundle| (self.selector)(bundle).is_empty())
    }

    /// Number of canonical phase entries, including phases whose selected
    /// property cache is currently empty.
    pub fn phase_count(&self) -> usize {
        self.phase_bundles.len()
    }

    /// Compatibility alias for `phase_count`.
    pub fn len(&self) -> usize {
        self.phase_count()
    }

    /// Number of phases that currently contain at least one value for this
    /// property. This differs from `phase_count` after cache invalidation.
    pub fn populated_phase_count(&self) -> usize {
        self.phase_bundles
            .values()
            .filter(|bundle| !(self.selector)(bundle).is_empty())
            .count()
    }

    /// Iterates without allocating or cloning phase identities. Prefer this
    /// for internal reads and hot-path adapters that do not need owned keys.
    pub fn iter_ref(&'a self) -> PhaseThermoPropertyIter<'a, T> {
        PhaseThermoPropertyIter {
            inner: self.phase_bundles.iter(),
            selector: self.selector,
        }
    }

    /// Compatibility iterator that yields owned phase keys. New code should
    /// prefer `iter_ref` when it can borrow the view.
    pub fn iter(
        &'a self,
    ) -> Box<dyn Iterator<Item = (Option<String>, &'a HashMap<String, T>)> + 'a> {
        Box::new(
            self.phase_bundles
                .iter()
                .map(|(phase, bundle)| (phase.clone(), (self.selector)(bundle))),
        )
    }

    pub fn to_owned_map(&self) -> HashMap<Option<String>, HashMap<String, T>>
    where
        T: Clone,
    {
        self.phase_bundles
            .iter()
            .map(|(phase, bundle)| (phase.clone(), (self.selector)(bundle).clone()))
            .collect()
    }
}

impl<'a, T> std::fmt::Debug for PhaseThermoPropertyView<'a, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PhaseThermoPropertyView")
            .field("phase_count", &self.phase_count())
            .field("populated_phase_count", &self.populated_phase_count())
            .finish()
    }
}

/// Borrowed view over nested phase caches without cloning the underlying data.
#[derive(Debug, Clone, Copy)]
pub enum NestedPhaseCacheView<'a, T> {
    Single(&'a HashMap<String, T>),
    Multi(PhaseThermoPropertyView<'a, T>),
}

static SINGLE_PHASE_KEY: Option<String> = None;

/// Allocation-free iterator over a nested cache view. The single-phase case
/// exposes one shared `None` phase key; the multi-phase case delegates to the
/// canonical borrowed phase iterator.
pub enum NestedPhaseCacheIter<'a, T> {
    Single(Option<&'a HashMap<String, T>>),
    Multi(PhaseThermoPropertyIter<'a, T>),
}

impl<'a, T: 'a> Iterator for NestedPhaseCacheIter<'a, T> {
    type Item = (&'a Option<String>, &'a HashMap<String, T>);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Single(value) => value.take().map(|map| (&SINGLE_PHASE_KEY, map)),
            Self::Multi(iter) => iter.next(),
        }
    }
}

impl<'a, T> NestedPhaseCacheView<'a, T> {
    pub fn get(&self, phase: &Option<String>) -> Option<&HashMap<String, T>> {
        match self {
            NestedPhaseCacheView::Single(map) => {
                if phase.is_none() {
                    Some(map)
                } else {
                    None
                }
            }
            NestedPhaseCacheView::Multi(map) => map.get(phase),
        }
    }

    pub fn iter(
        &'a self,
    ) -> Box<dyn Iterator<Item = (Option<String>, &'a HashMap<String, T>)> + 'a> {
        match self {
            NestedPhaseCacheView::Single(map) => Box::new(std::iter::once((None, *map))),
            NestedPhaseCacheView::Multi(map) => {
                Box::new(map.iter().map(|(phase, data)| (phase.clone(), data)))
            }
        }
    }

    /// Iterates phase keys and property maps by reference, avoiding phase-key
    /// cloning and iterator allocation. Prefer this for internal readers.
    pub fn iter_ref(&'a self) -> NestedPhaseCacheIter<'a, T> {
        match self {
            NestedPhaseCacheView::Single(map) => NestedPhaseCacheIter::Single(Some(*map)),
            NestedPhaseCacheView::Multi(map) => NestedPhaseCacheIter::Multi(map.iter_ref()),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            NestedPhaseCacheView::Single(map) => map.is_empty(),
            NestedPhaseCacheView::Multi(map) => map.is_empty(),
        }
    }

    /// Number of canonical phase entries represented by this view. A
    /// single-phase view always represents exactly one `None` phase, even when
    /// its selected property cache is empty.
    pub fn phase_count(&self) -> usize {
        match self {
            NestedPhaseCacheView::Single(_) => 1,
            NestedPhaseCacheView::Multi(map) => map.phase_count(),
        }
    }

    /// Number of phases with non-empty values for the selected property.
    pub fn populated_phase_count(&self) -> usize {
        match self {
            NestedPhaseCacheView::Single(map) => usize::from(!map.is_empty()),
            NestedPhaseCacheView::Multi(map) => map.populated_phase_count(),
        }
    }
}

/// Typed snapshot of a phase cache aligned with a layout revision.
#[derive(Debug, Clone, Copy)]
pub struct ThermoCacheSnapshot<'a, T> {
    pub layout_revision: usize,
    pub values: NestedPhaseCacheView<'a, T>,
}

impl<'a, T> ThermoCacheSnapshot<'a, T> {
    pub fn new(layout_revision: usize, values: NestedPhaseCacheView<'a, T>) -> Self {
        Self {
            layout_revision,
            values,
        }
    }
}

impl<'a> ThermoCacheSnapshot<'a, f64> {
    /// Projects a phase-keyed numeric cache into solver order without exposing
    /// its map iteration order. The snapshot and layout must originate from
    /// the same structural revision.
    pub fn ordered_values(
        &self,
        layout: &SystemLayout,
        layout_revision: usize,
        property_name: &str,
    ) -> SubsDataResult<Vec<f64>> {
        if self.layout_revision != layout_revision {
            return Err(SubsDataError::calculation_failed(
                property_name.to_string(),
                "ordered cache projection",
                format!(
                    "cache revision {} does not match layout revision {}",
                    self.layout_revision, layout_revision
                ),
            ));
        }

        layout
            .components()
            .iter()
            .map(|component| {
                let value = self
                    .values
                    .get(component.phase.as_option())
                    .and_then(|phase_values| phase_values.get(&component.substance))
                    .copied()
                    .ok_or_else(|| SubsDataError::MissingData {
                        field: format!("{property_name} cache"),
                        substance: component.label(),
                    })?;
                if !value.is_finite() {
                    return Err(SubsDataError::calculation_failed(
                        component.label(),
                        "ordered cache projection",
                        format!("{property_name} cache contains a non-finite value"),
                    ));
                }
                Ok(value)
            })
            .collect()
    }
}

pub(crate) fn build_cache_snapshot<'a, T>(
    layout_revision: usize,
    values: NestedPhaseCacheView<'a, T>,
) -> ThermoCacheSnapshot<'a, T> {
    ThermoCacheSnapshot::new(layout_revision, values)
}

#[inline]
pub(crate) fn build_single_cache_snapshot<'a, T>(
    layout_revision: usize,
    values: &'a HashMap<String, T>,
) -> ThermoCacheSnapshot<'a, T> {
    build_cache_snapshot(layout_revision, NestedPhaseCacheView::Single(values))
}

#[inline]
pub(crate) fn build_multi_cache_snapshot<'a, T>(
    layout_revision: usize,
    values: PhaseThermoPropertyView<'a, T>,
) -> ThermoCacheSnapshot<'a, T> {
    build_cache_snapshot(layout_revision, NestedPhaseCacheView::Multi(values))
}

/// Typed full-state snapshot for the phase system caches.
pub struct ThermoStateSnapshot<'a> {
    pub layout_revision: usize,
    pub dG: ThermoCacheSnapshot<'a, f64>,
    pub dG_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    pub dG_sym: ThermoCacheSnapshot<'a, Expr>,
    pub dS: ThermoCacheSnapshot<'a, f64>,
    pub dS_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    pub dS_sym: ThermoCacheSnapshot<'a, Expr>,
}

#[inline]
pub(crate) fn build_thermo_state_snapshot<'a>(
    layout_revision: usize,
    dG: ThermoCacheSnapshot<'a, f64>,
    dG_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    dG_sym: ThermoCacheSnapshot<'a, Expr>,
    dS: ThermoCacheSnapshot<'a, f64>,
    dS_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    dS_sym: ThermoCacheSnapshot<'a, Expr>,
) -> ThermoStateSnapshot<'a> {
    ThermoStateSnapshot {
        layout_revision,
        dG,
        dG_fun,
        dG_sym,
        dS,
        dS_fun,
        dS_sym,
    }
}

/// Typed thermodynamic result snapshot aligned to a canonical layout.
#[derive(Clone)]
pub struct ThermoResultSnapshot<'a> {
    pub layout_revision: usize,
    pub layout: SystemLayout,
    pub temperature: Option<f64>,
    pub pressure: Option<f64>,
    pub dG: ThermoCacheSnapshot<'a, f64>,
    pub dG_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    pub dG_sym: ThermoCacheSnapshot<'a, Expr>,
    pub dS: ThermoCacheSnapshot<'a, f64>,
    pub dS_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    pub dS_sym: ThermoCacheSnapshot<'a, Expr>,
}

impl<'a> ThermoResultSnapshot<'a> {
    pub fn component_labels(&self) -> Vec<String> {
        self.layout.component_labels()
    }

    pub fn phase_labels(&self) -> Vec<Option<String>> {
        self.layout
            .phases()
            .iter()
            .map(|phase| phase.as_option().clone())
            .collect()
    }

    pub fn component_count(&self) -> usize {
        self.layout.component_count()
    }

    pub fn phase_count(&self) -> usize {
        self.layout.phases().len()
    }

    /// Gibbs values in the canonical component order used by solver vectors.
    pub fn ordered_gibbs_values(&self) -> SubsDataResult<Vec<f64>> {
        self.dG
            .ordered_values(&self.layout, self.layout_revision, "Gibbs free energy")
    }

    /// Entropy values in the canonical component order used by solver vectors.
    pub fn ordered_entropy_values(&self) -> SubsDataResult<Vec<f64>> {
        self.dS
            .ordered_values(&self.layout, self.layout_revision, "entropy")
    }
}

#[inline]
pub(crate) fn build_thermo_result_snapshot<'a>(
    layout_revision: usize,
    layout: SystemLayout,
    temperature: Option<f64>,
    pressure: Option<f64>,
    dG: ThermoCacheSnapshot<'a, f64>,
    dG_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    dG_sym: ThermoCacheSnapshot<'a, Expr>,
    dS: ThermoCacheSnapshot<'a, f64>,
    dS_fun: ThermoCacheSnapshot<'a, PhaseFunction>,
    dS_sym: ThermoCacheSnapshot<'a, Expr>,
) -> ThermoResultSnapshot<'a> {
    ThermoResultSnapshot {
        layout_revision,
        layout,
        temperature,
        pressure,
        dG,
        dG_fun,
        dG_sym,
        dS,
        dS_fun,
        dS_sym,
    }
}

#[inline]
pub(crate) fn substitute_exprs<'a, I>(exprs: I, variable: &str, value: f64)
where
    I: IntoIterator<Item = &'a mut Expr>,
{
    for expr in exprs {
        *expr = expr.set_variable(variable, value).simplify();
    }
}
