//! # User Reactions Module - High-Level Kinetics API
//!
//! ## Purpose
//! This module provides the main user-facing API for the entire Kinetics module. It aggregates
//! functionality from all other kinetics modules into a single, comprehensive interface for
//! processing chemical reaction data. Acts as the primary entry point for kinetic analysis.
//!
//! ## Main Data Structures
//! - `KinData`: Central aggregation structure containing all reaction analysis results
//!   - `shortcut_reactions`: Vector of reaction shortcut names (e.g., "C1", "NUIG_42")
//!   - `reaction_ids`: Canonical reaction identities used as the internal source of truth
//!   - `vec_of_pairs`: Vector of (library, reaction_id) tuples for direct access
//!   - `vec_of_reaction_data`: Parsed reaction data as ReactionData structures
//!   - `vec_of_equations`: Human-readable reaction equation strings
//!   - `substances`: Vector of all unique substance names found
//!   - `stecheodata`: Complete stoichiometric analysis from StoichAnalyzer
//!
//! ## Key Logic Implementation
//! 1. **Reaction Collection**: Multiple pathways to gather reactions (shortcuts, ranges, mechanism construction)
//! 2. **Data Integration**: Combines kinetic parameters, stoichiometry, and molecular data
//! 3. **Shortcut Processing**: Converts user-friendly shortcuts ("C1..C10") to library addresses
//! 4. **Mechanism Construction**: Automatic reaction network generation from seed substances
//! 5. **Unified Analysis**: Single interface for all kinetic calculations and data export
//!
//! ## Usage Patterns
//! ```rust, ignore
//! // Method 1: Using reaction shortcuts
//! let mut kd = KinData::new();
//! kd.set_reactions_from_shortcut_range("C1..C10".to_string());
//! kd.get_reactions_from_shortcuts();
//! kd.kinetic_main();
//!
//! // Method 2: Mechanism construction
//! let mut kd = KinData::new();
//! kd.construct_mechanism(vec!["O".to_string(), "H2".to_string()], "NUIG".to_string());
//! kd.kinetic_main();
//!
//! // Method 3: Direct reaction input
//! let mut kd = KinData::new();
//! kd.set_reactions_directly(vec!["H + O2 = OH + O".to_string()], None);
//! kd.kinetic_main();
//! ```
//!
//! ## KinData Methods
//! ### Initialization
//! - `new()`: Create new empty KinData instance
//!
//! ### Reaction Input Methods
//! - `set_reactions_directly()`: Input reactions as equation strings with optional chemical groups
//! - `set_reactions_from_shortcut_range()`: Generate reaction shortcuts from range ("C1..C10")
//! - `get_reactions_from_shortcuts()`: Resolve shortcuts to actual reaction data from libraries
//! - `construct_mechanism()`: Auto-generate reaction network from seed substances and library
//!
//! ### Data Manipulation
//! - `append_reaction()`: Add reaction data as serde Values
//! - `append_reaction_with_shortcut()`: Add reactions with their shortcut names
//! - `remove_by_index()`: Remove reaction by index position
//! - `remove_reaction_by_eq()`: Remove reaction by equation string
//!
//! ### Analysis Methods
//! - `reactdata_parsing()`: Parse serde Values into ReactionData structures
//! - `analyze_reactions()`: Generate stoichiometric matrices and molecular data
//! - `kinetic_main()`: Combined parsing and analysis in one call
//!
//! ### Kinetic Calculations
//! - `calc_K_const_for_1_react()`: Calculate rate constant for single reaction
//! - `calc_K_const_for_all_reactions()`: Calculate and sort all rate constants
//! - `calc_K_const_for_all_reactions_forTrange()`: Calculate max rate constants over temperature range
//! -  `calc_sym_constants()`: Calculate symbolic constants for all reactions
//!
//! ### I/O Methods
//! - `save_raw_reactions()`: Export raw reaction data to JSON
//! - `save_reactions_with_shortcuts()`: Export with shortcut mappings
//! - `create_kinetics_document()`: Create structured kinetics document
//! - `load_reactions_from_json()`: Load reactions from JSON file
//! - `KinDataBuilder`: Builder-style API for constructing normalized KinData values
//!
//! ### Display Methods
//! - `pretty_print_substances_verbose()`: Print stoichiometric and elemental matrices
//! - `pretty_print_reaction_data()`: Print kinetic parameters in table format
//! - `pretty_print_kindata()`: Print both reaction data and matrices
//!
//! ## Interesting Features
//! - **Multi-Modal Input**: Supports shortcuts, ranges, direct equations, and automatic mechanism generation
//! - **Library Agnostic**: Works with any kinetic database (NUIG, Cantera, Aramco, etc.)
//! - **Integrated Analysis**: Combines kinetic parameters, stoichiometry, and molecular properties
//! - **Export Capabilities**: Built-in pretty printing and data serialization
//! - **Flexible Workflow**: Modular design allows partial analysis or complete processing
//! - **Error Recovery**: Robust handling of missing data and parsing errors

//mod mechfinder_api;
use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::Kinetics::kinetics_lib_api::KineticData;
use crate::Kinetics::mechfinder_api::{
    Mechanism_search, ReactionData, ReactionKinetics, parse_kinetic_data_vec,
};
use crate::Kinetics::parsetask::{
    decipher_vector_of_shortcuts, decipher_vector_of_shortcuts_to_pairs,
};
use crate::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use prettytable::{Cell, Row, Table, row};
use serde_json::json;
use serde_json::{Value, from_reader, to_writer_pretty};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};

use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;

use log::info;

///structure to store user task and reaction data
#[derive(Debug, Clone)]
pub struct KinData {
    /// canonical identity of each reaction; shortcuts remain optional aliases
    pub(crate) reaction_ids: Option<Vec<ReactionIdentity>>,
    /// vector of reaction shortcut names
    pub(crate) shortcut_reactions: Option<Vec<String>>,
    /// Compatibility cache for paired library/reaction identifiers.
    /// The canonical read path now derives these pairs from reaction identities.
    pub(crate) vec_of_pairs: Option<Vec<(String, String)>>,
    /// vector of reaction data in the form of serde Values
    pub(crate) vec_of_reaction_Values: Option<Vec<Value>>,
    /// data of all reactions
    pub(crate) vec_of_reaction_data: Option<Vec<ReactionData>>,
    /// vector of equations of reactions
    pub(crate) vec_of_equations: Vec<String>,
    /// vector of substance names
    pub(crate) substances: Vec<String>,
    /// Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
    pub(crate) groups: Option<HashMap<String, HashMap<String, usize>>>,
    /// matrix of stoichiometric coefficients and other matrices
    pub(crate) stecheodata: StoichAnalyzer,
    /// vector of symbolic reaction constants
    pub(crate) K_sym_vec: Option<Vec<Expr>>,
    /// lifecycle state of the reaction set
    pub(crate) state: KinDataState,
    /// Normalized per-reaction view used by export, formatting, and sorted workflows.
    pub(crate) every_reaction: Option<Vec<EveryReaction>>,
}
#[derive(Debug, Clone)]
pub struct EveryReaction {
    pub reaction_id: ReactionIdentity,
    pub shortcut: Option<String>,
    pub lib_and_id: Option<(String, String)>,
    pub reaction: ReactionData,
    pub equation: String,
    pub K_sym: Option<Expr>,
}
#[derive(Debug, Clone)]
struct ReactionBranchSnapshot {
    reaction_ids: Vec<ReactionIdentity>,
    shortcut_reactions: Option<Vec<String>>,
    vec_of_pairs: Option<Vec<(String, String)>>,
    vec_of_reaction_values: Option<Vec<Value>>,
    vec_of_reaction_data: Option<Vec<ReactionData>>,
    vec_of_equations: Vec<String>,
    substances: Vec<String>,
    rebuild_every_reaction: bool,
}

#[derive(Debug, Clone)]
struct ReorderedReactionViews {
    equations: Vec<String>,
    shortcuts: Option<Vec<String>>,
    pairs: Option<Vec<(String, String)>>,
    reaction_values: Option<Vec<Value>>,
    reaction_data: Vec<ReactionData>,
    k_sym: Option<Vec<Expr>>,
    reaction_ids: Option<Vec<ReactionIdentity>>,
}

#[derive(Debug, Clone)]
struct NormalizedReactionSnapshot {
    reaction_data: Vec<ReactionData>,
    reaction_ids: Vec<ReactionIdentity>,
    equations: Vec<String>,
    reaction_values: Option<Vec<Value>>,
    k_sym: Option<Vec<Expr>>,
}

impl NormalizedReactionSnapshot {
    /// Reorder a normalized reaction snapshot while keeping every aligned view together.
    fn reordered(&self, indices: &[usize]) -> KineticsResult<ReorderedReactionViews> {
        if indices.len() != self.reaction_data.len() {
            return Err(KineticsError::LengthMismatch {
                context: "KinData::build_reordered_reaction_views(indices vs vec_of_reaction_data)",
                left: indices.len(),
                right: self.reaction_data.len(),
            });
        }
        if indices
            .iter()
            .any(|&index| index >= self.reaction_data.len())
        {
            return Err(KineticsError::InvalidReactionData(
                "reorder indices point outside the reaction data range".to_string(),
            ));
        }

        let equations = indices
            .iter()
            .map(|&i| self.equations[i].clone())
            .collect::<Vec<_>>();
        let reaction_data = indices
            .iter()
            .map(|&i| self.reaction_data[i].clone())
            .collect::<Vec<_>>();
        let reaction_ids = indices
            .iter()
            .map(|&i| self.reaction_ids[i].clone())
            .collect::<Vec<_>>();

        let mut reaction_values = self.reaction_values.clone();
        KinData::reorder_option_vec(&mut reaction_values, indices);

        let mut k_sym = self.k_sym.clone();
        KinData::reorder_option_vec(&mut k_sym, indices);

        let shortcuts = KinData::derived_shortcut_names_from_ids(&reaction_ids);
        let pairs = KinData::derived_library_pairs_from_ids(&reaction_ids);

        Ok(ReorderedReactionViews {
            equations,
            shortcuts,
            pairs,
            reaction_values,
            reaction_data,
            k_sym,
            reaction_ids: Some(reaction_ids),
        })
    }
}

impl ReactionBranchSnapshot {
    /// Build a direct-equation branch from raw equation strings.
    fn direct_equations(reactions: Vec<String>) -> Self {
        let reaction_ids = reactions
            .iter()
            .enumerate()
            .map(|(index, _)| ReactionIdentity::Document {
                source: DocumentSourceKind::Direct.as_str().to_string(),
                index,
            })
            .collect();
        Self {
            reaction_ids,
            shortcut_reactions: None,
            vec_of_pairs: None,
            vec_of_reaction_values: None,
            vec_of_reaction_data: None,
            vec_of_equations: reactions,
            substances: Vec::new(),
            rebuild_every_reaction: false,
        }
    }

    /// Build a shortcut-selection branch without hydrated reaction data.
    fn shortcut_selection(shortcuts: Vec<String>) -> Self {
        Self {
            reaction_ids: shortcuts
                .iter()
                .cloned()
                .map(ReactionIdentity::Shortcut)
                .collect(),
            shortcut_reactions: Some(shortcuts),
            vec_of_pairs: None,
            vec_of_reaction_values: None,
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            rebuild_every_reaction: false,
        }
    }

    /// Build a shortcut-loaded branch with raw reaction values and shortcut aliases.
    fn shortcut_values(shortcuts: Vec<String>, reaction_values: Vec<Value>) -> Self {
        Self {
            reaction_ids: shortcuts
                .iter()
                .cloned()
                .map(ReactionIdentity::Shortcut)
                .collect(),
            shortcut_reactions: Some(shortcuts),
            vec_of_pairs: None,
            vec_of_reaction_values: Some(reaction_values),
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            rebuild_every_reaction: false,
        }
    }

    /// Build a raw append branch with document-style reaction ids.
    fn raw_values(
        reaction_values: Vec<Value>,
        source: DocumentSourceKind,
        start_index: usize,
    ) -> Self {
        let reaction_ids = (0..reaction_values.len())
            .map(|offset| ReactionIdentity::Document {
                source: source.as_str().to_string(),
                index: start_index + offset,
            })
            .collect();
        Self {
            reaction_ids,
            shortcut_reactions: None,
            vec_of_pairs: None,
            vec_of_reaction_values: Some(reaction_values),
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            rebuild_every_reaction: false,
        }
    }

    /// Build a parsed-data branch when callers already have typed ReactionData rows.
    fn parsed_reaction_data(vec_of_reaction_data: Vec<ReactionData>) -> Self {
        let vec_of_equations = vec_of_reaction_data
            .iter()
            .map(|reaction| reaction.eq.clone())
            .collect::<Vec<_>>();
        let reaction_ids = vec_of_reaction_data
            .iter()
            .enumerate()
            .map(|(index, _)| ReactionIdentity::Document {
                source: DocumentSourceKind::Direct.as_str().to_string(),
                index,
            })
            .collect();

        Self {
            reaction_ids,
            shortcut_reactions: None,
            vec_of_pairs: None,
            vec_of_reaction_values: None,
            vec_of_reaction_data: Some(vec_of_reaction_data),
            vec_of_equations,
            substances: Vec::new(),
            rebuild_every_reaction: true,
        }
    }

    /// Build a library-resolved branch with reaction values but without parsed equations.
    fn library_resolved(
        shortcut_reactions: Option<Vec<String>>,
        vec_of_pairs: Vec<(String, String)>,
        vec_of_reaction_values: Option<Vec<Value>>,
    ) -> Self {
        Self {
            reaction_ids: vec_of_pairs
                .iter()
                .map(|(lib, reaction_id)| ReactionIdentity::LibraryReaction {
                    library: lib.clone(),
                    reaction_id: reaction_id.clone(),
                })
                .collect(),
            shortcut_reactions,
            vec_of_pairs: Some(vec_of_pairs),
            vec_of_reaction_values,
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            rebuild_every_reaction: false,
        }
    }

    /// Build a mechanism branch from the generated mechanism reaction rows.
    fn mechanism_built(
        task_library: String,
        reactions: Vec<String>,
        vec_of_reaction_data: Vec<ReactionData>,
        vec_of_equations: Vec<String>,
        substances: Vec<String>,
    ) -> Self {
        let shortcut_reactions = reactions
            .iter()
            .map(|reaction| format!("{}_{}", task_library, reaction))
            .collect();
        let vec_of_pairs: Vec<(String, String)> = reactions
            .iter()
            .map(|reaction| (task_library.clone(), reaction.clone()))
            .collect();
        Self {
            reaction_ids: reactions
                .iter()
                .map(|reaction| ReactionIdentity::LibraryReaction {
                    library: task_library.clone(),
                    reaction_id: reaction.clone(),
                })
                .collect(),
            shortcut_reactions: Some(shortcut_reactions),
            vec_of_pairs: Some(vec_of_pairs),
            vec_of_reaction_values: None,
            vec_of_reaction_data: Some(vec_of_reaction_data),
            vec_of_equations,
            substances,
            rebuild_every_reaction: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ReactionIdentity {
    Shortcut(String),
    LibraryReaction {
        library: String,
        reaction_id: String,
    },
    Document {
        source: String,
        index: usize,
    },
}

impl ReactionIdentity {
    /// Return the internal document source tag when this identity refers to a document branch.
    fn document_source_kind(&self) -> Option<DocumentSourceKind> {
        let ReactionIdentity::Document { source, .. } = self else {
            return None;
        };

        match source.as_str() {
            "direct" => Some(DocumentSourceKind::Direct),
            "raw" => Some(DocumentSourceKind::Raw),
            "append" => Some(DocumentSourceKind::Append),
            _ => None,
        }
    }

    /// Report whether this identity belongs to the shortcut branch.
    fn is_shortcut(&self) -> bool {
        matches!(self, ReactionIdentity::Shortcut(_))
    }

    /// Report whether this identity belongs to the library-backed branch.
    fn is_library_reaction(&self) -> bool {
        matches!(self, ReactionIdentity::LibraryReaction { .. })
    }

    /// Report whether this identity belongs to a document/raw branch.
    fn is_document_branch(&self) -> bool {
        self.document_source_kind().is_some()
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum KinDataState {
    #[default]
    Empty,
    ShortcutsSelected,
    LibraryResolved,
    DirectReactionsLoaded,
    ReactionDataParsed,
    MechanismBuilt,
    Analyzed,
    Sorted,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReactionBranchKind {
    Empty,
    Shortcut,
    Library,
    Direct,
    Parsed,
    Mechanism,
    Normalized,
}

/// Internal snapshot of the current workflow shape.
#[derive(Debug, Clone, Copy)]
struct WorkflowSnapshot {
    branch_kind: ReactionBranchKind,
    normalized_view_ready: bool,
    has_parsed_branch: bool,
    has_equations: bool,
    has_mechanism_branch: bool,
    has_analyzed_artifacts: bool,
}

impl WorkflowSnapshot {
    /// Convert the current workflow shape into the most specific visible state.
    fn visible_state(self) -> KinDataState {
        if self.has_mechanism_branch {
            KinDataState::MechanismBuilt
        } else if self.normalized_view_ready {
            KinDataState::Sorted
        } else if self.has_analyzed_artifacts {
            KinDataState::Analyzed
        } else {
            KinData::state_from_branch_kind(self.branch_kind)
        }
    }
}

/// Internal document source tag used for raw and direct reaction payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DocumentSourceKind {
    Direct,
    Raw,
    Append,
}

impl DocumentSourceKind {
    /// Convert the internal tag to the persisted string stored in `ReactionIdentity`.
    fn as_str(self) -> &'static str {
        match self {
            DocumentSourceKind::Direct => "direct",
            DocumentSourceKind::Raw => "raw",
            DocumentSourceKind::Append => "append",
        }
    }

    /// Report whether a raw append can extend the current document branch.
    fn can_extend(self) -> bool {
        matches!(self, DocumentSourceKind::Append)
    }
}
#[derive(Debug, Clone)]
pub struct KinDataBuilder {
    inner: KinData,
}
impl Default for KinDataBuilder {
    fn default() -> Self {
        Self::new()
    }
}
impl KinDataBuilder {
    pub fn new() -> Self {
        Self {
            inner: KinData::new(),
        }
    }

    fn reset_reaction_input(&mut self) {
        // Builder setters should not leak stale state from a previous workflow branch.
        self.inner.reset_to_empty();
    }

    /// Add a direct-reaction branch to the builder.
    pub fn with_direct_reactions(mut self, reactions: Vec<String>) -> KineticsResult<Self> {
        // Start from raw equations and let the inner API derive the canonical ids.
        self.reset_reaction_input();
        self.inner.set_reactions_directly(reactions, None)?;
        Ok(self)
    }

    /// Add a shortcut branch to the builder.
    pub fn with_shortcuts(mut self, shortcuts: Vec<String>) -> KineticsResult<Self> {
        // Keep shortcuts as the user-facing input, but store them through canonical ids too.
        self.reset_reaction_input();
        self.inner
            .install_reaction_branch(ReactionBranchSnapshot::shortcut_selection(shortcuts))?;
        Ok(self)
    }

    /// Expand and add a shortcut range to the builder.
    pub fn with_shortcut_range(mut self, shortcut_range: String) -> KineticsResult<Self> {
        // Parse the range once, then reuse the same shortcut normalization path.
        self.reset_reaction_input();
        self.inner
            .set_reactions_from_shortcut_range(shortcut_range)?;
        Ok(self)
    }

    /// Add a library-backed branch to the builder.
    pub fn with_library_reactions(
        mut self,
        library: String,
        ids: Vec<String>,
    ) -> KineticsResult<Self> {
        // Library workflows must start from a clean branch so stale direct or shortcut data do not leak in.
        self.reset_reaction_input();
        // Library pairs are stored both as canonical ids and as the export-friendly pair view.
        let pairs: Vec<(String, String)> = ids
            .into_iter()
            .map(|reaction_id| (library.clone(), reaction_id))
            .collect();
        self.inner
            .install_reaction_branch(ReactionBranchSnapshot::library_resolved(None, pairs, None))?;
        Ok(self)
    }

    /// Attach optional atom-group aliases to the builder.
    pub fn with_groups(mut self, groups: Option<HashMap<String, HashMap<String, usize>>>) -> Self {
        // Preserve custom group aliases alongside the chosen reaction source.
        self.inner.groups = groups;
        self
    }

    /// Backward-compatible alias for `with_direct_reactions`.
    pub fn direct_reactions(mut self, reactions: Vec<String>) -> KineticsResult<Self> {
        self = self.with_direct_reactions(reactions)?;
        Ok(self)
    }

    /// Backward-compatible alias for `with_shortcuts`.
    pub fn shortcut_reactions(mut self, shortcuts: Vec<String>) -> KineticsResult<Self> {
        self = self.with_shortcuts(shortcuts)?;
        Ok(self)
    }

    /// Backward-compatible alias for `with_shortcut_range`.
    pub fn shortcut_range(mut self, shortcut_range: String) -> KineticsResult<Self> {
        self = self.with_shortcut_range(shortcut_range)?;
        Ok(self)
    }

    /// Backward-compatible alias for `with_library_reactions`.
    pub fn library_reactions(mut self, library: String, ids: Vec<String>) -> KineticsResult<Self> {
        self = self.with_library_reactions(library, ids)?;
        Ok(self)
    }

    /// Backward-compatible alias for `with_groups`.
    pub fn groups(mut self, groups: Option<HashMap<String, HashMap<String, usize>>>) -> Self {
        self = self.with_groups(groups);
        self
    }

    pub fn build(mut self) -> KineticsResult<KinData> {
        // Normalize the assembled data before returning it to the caller.
        self.inner.normalize_before_use()?;
        Ok(self.inner)
    }
}
impl KinData {
    pub fn builder() -> KinDataBuilder {
        KinDataBuilder::new()
    }

    /// Build a direct-reaction branch in one call.
    pub fn from_direct_reactions(reactions: Vec<String>) -> KineticsResult<Self> {
        Self::builder().with_direct_reactions(reactions)?.build()
    }

    /// Build a shortcut branch in one call.
    pub fn from_shortcuts(shortcuts: Vec<String>) -> KineticsResult<Self> {
        Self::builder().with_shortcuts(shortcuts)?.build()
    }

    /// Build a shortcut-range branch in one call.
    pub fn from_shortcut_range(shortcut_range: String) -> KineticsResult<Self> {
        Self::builder().with_shortcut_range(shortcut_range)?.build()
    }

    /// Build a library-backed branch in one call.
    pub fn from_library_reactions(library: String, ids: Vec<String>) -> KineticsResult<Self> {
        Self::builder()
            .with_library_reactions(library, ids)?
            .build()
    }

    /// Build a branch and attach optional group aliases in one call.
    pub fn from_direct_reactions_with_groups(
        reactions: Vec<String>,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> KineticsResult<Self> {
        Self::builder()
            .with_direct_reactions(reactions)?
            .with_groups(groups)
            .build()
    }

    /// Return the current workflow state.
    pub fn workflow_state(&self) -> KinDataState {
        self.state
    }

    /// Return the current workflow state using a shorter, more idiomatic name.
    pub fn state(&self) -> KinDataState {
        self.state
    }

    /// Return the canonical reaction identities, if available.
    pub fn canonical_reaction_ids(&self) -> Option<&[ReactionIdentity]> {
        self.reaction_ids.as_deref()
    }

    /// Return the current shortcut aliases, if any.
    pub fn shortcut_names(&self) -> Option<Vec<String>> {
        self.derived_shortcut_names()
    }

    /// Return the library reaction pairs, if any.
    pub fn library_pairs(&self) -> Option<Vec<(String, String)>> {
        self.derived_library_pairs()
    }

    /// Build the reaction map grouped by library from the canonical reaction ids.
    pub fn reaction_map(&self) -> Option<HashMap<String, Vec<String>>> {
        self.derived_reaction_map()
    }

    /// Return the raw reaction payload values, if any.
    pub fn reaction_values(&self) -> Option<&[Value]> {
        self.vec_of_reaction_Values.as_deref()
    }

    /// Report whether raw reaction payload values are available.
    pub fn has_reaction_values(&self) -> bool {
        self.reaction_values()
            .is_some_and(|values| !values.is_empty())
    }

    /// Return the parsed reaction data, if any.
    pub fn reaction_data(&self) -> Option<&[ReactionData]> {
        self.vec_of_reaction_data.as_deref()
    }

    /// Report whether parsed reaction data is available.
    pub fn has_reaction_data(&self) -> bool {
        self.reaction_data().is_some_and(|data| !data.is_empty())
    }

    /// Return the human-readable reaction equations.
    pub fn equations(&self) -> &[String] {
        &self.vec_of_equations
    }

    /// Return the substance names produced by stoichiometric analysis.
    pub fn substances(&self) -> &[String] {
        &self.substances
    }

    /// Return optional custom group definitions used by formula parsing.
    pub fn groups(&self) -> Option<&HashMap<String, HashMap<String, usize>>> {
        self.groups.as_ref()
    }

    /// Return the symbolic constants, if any.
    pub fn symbolic_constants(&self) -> Option<&[Expr]> {
        self.K_sym_vec.as_deref()
    }

    /// Report whether symbolic constants are available.
    pub fn has_symbolic_constants(&self) -> bool {
        self.symbolic_constants()
            .is_some_and(|constants| !constants.is_empty())
    }

    /// Return the normalized reaction records, if any.
    pub fn normalized_reactions(&self) -> Option<&[EveryReaction]> {
        self.every_reaction.as_deref()
    }

    /// Return the stoichiometry analyzer backing the current facade.
    pub fn stoichiometric_analyzer(&self) -> &StoichAnalyzer {
        &self.stecheodata
    }

    /// Return the number of reactions currently represented by the facade.
    pub fn reaction_count(&self) -> usize {
        self.canonical_reaction_count().unwrap_or(0)
    }

    /// Report whether shortcut aliases are available.
    pub fn has_shortcuts(&self) -> bool {
        self.shortcut_names()
            .is_some_and(|shortcuts| !shortcuts.is_empty())
    }

    /// Report whether library pairs are available.
    pub fn has_library_pairs(&self) -> bool {
        self.library_pairs().is_some_and(|pairs| !pairs.is_empty())
    }

    /// Return whether the current branch is fully normalized.
    pub fn has_normalized_view(&self) -> bool {
        self.normalized_reactions()
            .is_some_and(|rows| !rows.is_empty())
    }

    /// Report whether the current workflow is in mechanism-built parsed form.
    pub fn is_mechanism_built(&self) -> bool {
        self.state == KinDataState::MechanismBuilt && self.has_reaction_data()
    }

    /// Report whether the facade is currently in a consistent normalized state.
    pub fn is_normalized(&self) -> bool {
        self.state == KinDataState::Sorted && self.has_normalized_view()
    }

    /// Report whether the facade is currently sorted.
    pub fn is_sorted(&self) -> bool {
        self.state == KinDataState::Sorted
    }

    /// Report whether the facade has parsed reaction data.
    pub fn is_parsed(&self) -> bool {
        matches!(
            self.state,
            KinDataState::ReactionDataParsed | KinDataState::Analyzed | KinDataState::Sorted
        ) && self.has_reaction_data()
    }

    /// Report whether the facade has already been analyzed.
    pub fn is_analyzed(&self) -> bool {
        matches!(self.state, KinDataState::Analyzed | KinDataState::Sorted)
            && self.has_analysis_results()
    }

    pub fn new() -> Self {
        Self {
            reaction_ids: None,
            shortcut_reactions: None,
            vec_of_pairs: None,
            vec_of_reaction_Values: None,
            vec_of_reaction_data: None,
            vec_of_equations: Vec::new(),
            substances: Vec::new(),
            groups: None,
            stecheodata: StoichAnalyzer::new(),
            K_sym_vec: None,
            state: KinDataState::Empty,
            every_reaction: None,
        }
    }
    fn clear_reaction_cache(&mut self) {
        // These fields are derived views and must be rebuilt after state changes.
        self.vec_of_reaction_data = None;
        self.vec_of_equations.clear();
        self.K_sym_vec = None;
        self.every_reaction = None;
    }

    /// Clear the branch metadata that identifies a reaction set but does not describe it.
    fn clear_reaction_branch_metadata(&mut self) {
        self.reaction_ids = None;
        self.shortcut_reactions = None;
        self.vec_of_pairs = None;
        self.vec_of_reaction_Values = None;
    }

    /// Clear the reaction-domain workspaces that are rebuilt from equations or parsed data.
    fn clear_reaction_domain_views(&mut self) {
        self.substances.clear();
        self.stecheodata = StoichAnalyzer::new();
    }

    fn reset_to_empty(&mut self) {
        // Reset the whole reaction workspace before starting a new user workflow.
        self.clear_reaction_branch_metadata();
        self.clear_reaction_cache();
        self.clear_reaction_domain_views();
        self.groups = None;
        self.transition_to_state(KinDataState::Empty);
    }

    /// Clear the current reaction branch while preserving external user settings like groups.
    fn reset_reaction_branch(&mut self) {
        // Branch switches must not keep any of the previous reaction payload or derived views.
        self.clear_reaction_branch_metadata();
        self.clear_reaction_cache();
        self.clear_reaction_domain_views();
        self.transition_to_state(KinDataState::Empty);
    }

    fn reaction_identity_from_metadata(&self, index: usize) -> ReactionIdentity {
        if let Some(ids) = self.canonical_reaction_ids() {
            if let Some(id) = ids.get(index) {
                // Canonical ids win whenever they are already available.
                return id.clone();
            }
        }
        if let Some(every_reaction) = self.normalized_reactions() {
            if let Some(reaction) = every_reaction.get(index) {
                // A materialized normalized result is richer than compatibility caches.
                return reaction.reaction_id.clone();
            }
        }
        // Legacy metadata is only a fallback when the canonical id cache has not been rebuilt yet.
        if let Some(pairs) = self.vec_of_pairs.as_ref() {
            if let Some((library, reaction_id)) = pairs.get(index) {
                return ReactionIdentity::LibraryReaction {
                    library: library.clone(),
                    reaction_id: reaction_id.clone(),
                };
            }
        }
        if let Some(shortcuts) = self.shortcut_reactions.as_ref() {
            if let Some(shortcut) = shortcuts.get(index) {
                return ReactionIdentity::Shortcut(shortcut.clone());
            }
        }
        if self.reaction_values().is_some() {
            return ReactionIdentity::Document {
                source: DocumentSourceKind::Raw.as_str().to_string(),
                index,
            };
        }
        ReactionIdentity::Document {
            source: DocumentSourceKind::Raw.as_str().to_string(),
            index,
        }
    }

    fn reaction_id_to_document_key(reaction_id: &ReactionIdentity) -> (String, String) {
        // Exports need a stable string key for every canonical reaction identity.
        match reaction_id {
            ReactionIdentity::Shortcut(shortcut) => ("SHORTCUTS".to_string(), shortcut.clone()),
            ReactionIdentity::LibraryReaction {
                library,
                reaction_id,
            } => (library.clone(), reaction_id.clone()),
            ReactionIdentity::Document { source, index } => (source.clone(), index.to_string()),
        }
    }

    /// Overwrite the live branch payload in one place before validation and state updates.
    fn assign_reaction_branch_fields(&mut self, branch: ReactionBranchSnapshot) {
        self.reaction_ids = Some(branch.reaction_ids);
        self.shortcut_reactions = branch.shortcut_reactions;
        self.vec_of_pairs = branch.vec_of_pairs;
        self.vec_of_reaction_Values = branch.vec_of_reaction_values;
        self.vec_of_reaction_data = branch.vec_of_reaction_data;
        self.vec_of_equations = branch.vec_of_equations;
        self.substances = branch.substances;
    }

    /// Overwrite the reordered runtime views in one place so sorting stays coherent.
    fn assign_reordered_reaction_views(&mut self, views: ReorderedReactionViews) {
        self.vec_of_equations = views.equations;
        self.shortcut_reactions = views.shortcuts;
        self.vec_of_pairs = views.pairs;
        self.vec_of_reaction_Values = views.reaction_values;
        self.vec_of_reaction_data = Some(views.reaction_data);
        self.K_sym_vec = views.k_sym;
        self.reaction_ids = views.reaction_ids;
    }

    /// Store symbolic constants in one place before the usual sync and state refresh.
    fn store_symbolic_constants(&mut self, vec_of_k_sym: Vec<Expr>) {
        self.K_sym_vec = Some(vec_of_k_sym);
    }

    /// Install a complete reaction branch in one place so branch switches remain transactional.
    fn install_reaction_branch(&mut self, branch: ReactionBranchSnapshot) -> KineticsResult<()> {
        let rebuild_every_reaction = branch.rebuild_every_reaction;
        self.reset_reaction_branch();
        self.assign_reaction_branch_fields(branch);
        self.finalize_reaction_mutation(rebuild_every_reaction, None)?;
        Ok(())
    }

    /// Write the lifecycle state after it has been validated against the current payload.
    fn transition_to_state(&mut self, state: KinDataState) {
        self.state = state;
    }

    /// Validate a workflow state and commit it only when the current payload is compatible.
    fn set_state_checked(&mut self, state: KinDataState) -> KineticsResult<()> {
        self.validate_state_contract(state)?;
        self.transition_to_state(state);
        Ok(())
    }

    /// Verify that the current workflow state is allowed for the requested operation.
    fn ensure_state_allowed(
        &self,
        operation: &'static str,
        allowed_states: &[KinDataState],
    ) -> KineticsResult<()> {
        if allowed_states.iter().any(|state| *state == self.state) {
            Ok(())
        } else {
            Err(KineticsError::InvalidState(format!(
                "{} is not allowed in state {:?}",
                operation, self.state
            )))
        }
    }

    /// Replace parsed reaction data and equation cache together.
    fn store_parsed_reaction_views(
        &mut self,
        vec_of_reaction_data: Vec<ReactionData>,
        vec_of_equations: Vec<String>,
    ) {
        self.vec_of_reaction_data = Some(vec_of_reaction_data);
        self.vec_of_equations = vec_of_equations;
    }

    /// Refresh the lightweight substance list from equations without calculating molar masses.
    ///
    /// This is intentionally weaker than full stoichiometric analysis. It is used by typed
    /// reaction-data inputs so downstream reactor code can validate boundary conditions even
    /// when formula-level element data will be supplied or checked later.
    fn refresh_substances_from_equations(&mut self) -> KineticsResult<()> {
        if self.equations().is_empty() {
            self.substances.clear();
            return Ok(());
        }

        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = self.equations().to_vec();
        analyzer.groups = self.groups.clone();
        analyzer.search_substances()?;
        self.substances = analyzer.substances;
        Ok(())
    }

    /// Report whether the materialized per-reaction result is ready for export/formatting.
    fn normalized_reaction_view_is_ready(&self) -> bool {
        let Some(every_reaction) = self.normalized_reactions() else {
            return false;
        };
        if every_reaction.is_empty() {
            return false;
        }

        // `reaction_ids` is a request/identity layer. `every_reaction` is the materialized
        // result layer. If both are present they must agree, but the result layer can stand
        // alone for read-only export paths.
        self.canonical_reaction_ids().map_or(true, |ids| {
            ids.len() == every_reaction.len()
                && ids
                    .iter()
                    .zip(every_reaction.iter())
                    .all(|(id, reaction)| id == &reaction.reaction_id)
        })
    }

    /// Detect shortcut-only workflows without peeking at naming conventions.
    fn has_shortcut_branch(&self) -> bool {
        self.canonical_reaction_ids()
            .is_some_and(|ids| !ids.is_empty() && ids.iter().all(ReactionIdentity::is_shortcut))
    }

    /// Detect library-backed workflows by their structural metadata.
    fn has_library_branch(&self) -> bool {
        self.canonical_reaction_ids().is_some_and(|ids| {
            !ids.is_empty() && ids.iter().all(ReactionIdentity::is_library_reaction)
        })
    }

    /// Detect direct-reaction workflows from the raw value payload or document ids.
    fn has_direct_branch(&self) -> bool {
        self.reaction_values()
            .is_some_and(|values| !values.is_empty())
            || self.canonical_reaction_ids().is_some_and(|ids| {
                !ids.is_empty() && ids.iter().all(ReactionIdentity::is_document_branch)
            })
    }

    /// Detect parsed reaction data that is ready for stoichiometric analysis.
    fn has_parsed_branch(&self) -> bool {
        self.reaction_data().is_some_and(|data| !data.is_empty())
    }

    /// Detect mechanism-built data as a parsed branch with richer source metadata.
    fn has_mechanism_branch(&self) -> bool {
        self.is_mechanism_built()
    }

    /// Report whether the stoichiometric analysis cache is populated and internally consistent.
    fn has_analysis_results(&self) -> bool {
        !self.stecheodata.reactions.is_empty() && !self.stecheodata.stecheo_matrx.is_empty()
    }

    /// Report whether any analyzed artifact is still present in the workspace.
    fn has_analyzed_artifacts(&self) -> bool {
        self.has_analysis_results()
            || self
                .symbolic_constants()
                .is_some_and(|k_sym| !k_sym.is_empty())
    }

    /// Derive the grouped reaction map from canonical ids only.
    fn derived_reaction_map(&self) -> Option<HashMap<String, Vec<String>>> {
        self.canonical_reaction_ids()
            .and_then(Self::derived_reaction_map_from_ids)
    }

    /// Derive the library/reaction pair view from canonical ids only.
    fn derived_library_pairs(&self) -> Option<Vec<(String, String)>> {
        self.canonical_reaction_ids()
            .and_then(Self::derived_library_pairs_from_ids)
    }

    /// Derive shortcut aliases from canonical ids only.
    fn derived_shortcut_names(&self) -> Option<Vec<String>> {
        self.canonical_reaction_ids()
            .and_then(Self::derived_shortcut_names_from_ids)
    }

    /// Derive grouped reaction metadata from canonical ids.
    fn derived_reaction_map_from_ids(
        ids: &[ReactionIdentity],
    ) -> Option<HashMap<String, Vec<String>>> {
        if ids.is_empty() {
            return None;
        }

        // Mixed canonical ids are a stale-state smell: each compatibility view should describe
        // one coherent branch, not a partial projection of several branches at once.
        if Self::ids_are_all_shortcuts(ids) {
            let shortcuts = Self::derived_shortcut_names_from_ids(ids)?;
            let shortcut_refs: Vec<&str> =
                shortcuts.iter().map(|shortcut| shortcut.as_str()).collect();
            return Some(decipher_vector_of_shortcuts(shortcut_refs));
        }

        if !Self::ids_are_all_library_reactions(ids) {
            return None;
        }

        let mut grouped: HashMap<String, Vec<String>> = HashMap::new();
        for reaction_id in ids {
            if let ReactionIdentity::LibraryReaction {
                library,
                reaction_id,
            } = reaction_id
            {
                grouped
                    .entry(library.clone())
                    .or_default()
                    .push(reaction_id.clone());
            }
        }

        if grouped.is_empty() {
            None
        } else {
            Some(grouped)
        }
    }

    /// Derive library/reaction pairs from canonical ids.
    fn derived_library_pairs_from_ids(ids: &[ReactionIdentity]) -> Option<Vec<(String, String)>> {
        if !Self::ids_are_all_library_reactions(ids) {
            return None;
        }

        let pairs: Vec<(String, String)> = ids
            .iter()
            .map(|reaction_id| match reaction_id {
                ReactionIdentity::LibraryReaction {
                    library,
                    reaction_id,
                } => (library.clone(), reaction_id.clone()),
                _ => unreachable!("library-only canonical view was validated beforehand"),
            })
            .collect();
        if pairs.is_empty() { None } else { Some(pairs) }
    }

    /// Derive shortcut aliases from canonical ids.
    fn derived_shortcut_names_from_ids(ids: &[ReactionIdentity]) -> Option<Vec<String>> {
        if !Self::ids_are_all_shortcuts(ids) {
            return None;
        }

        let shortcuts: Vec<String> = ids
            .iter()
            .map(|reaction_id| match reaction_id {
                ReactionIdentity::Shortcut(shortcut) => shortcut.clone(),
                _ => unreachable!("shortcut-only canonical view was validated beforehand"),
            })
            .collect();
        if shortcuts.is_empty() {
            None
        } else {
            Some(shortcuts)
        }
    }

    /// Check whether every canonical id belongs to the shortcut branch.
    fn ids_are_all_shortcuts(ids: &[ReactionIdentity]) -> bool {
        !ids.is_empty()
            && ids
                .iter()
                .all(|id| matches!(id, ReactionIdentity::Shortcut(_)))
    }

    /// Check whether every canonical id belongs to the library branch.
    fn ids_are_all_library_reactions(ids: &[ReactionIdentity]) -> bool {
        !ids.is_empty()
            && ids
                .iter()
                .all(|id| matches!(id, ReactionIdentity::LibraryReaction { .. }))
    }

    fn remove_option_index<T>(slot: &mut Option<Vec<T>>, index: usize) {
        if let Some(values) = slot.as_mut() {
            if index < values.len() {
                values.remove(index);
            }
        }
    }

    fn reorder_option_vec<T: Clone>(slot: &mut Option<Vec<T>>, indices: &[usize]) {
        if let Some(values) = slot.as_mut() {
            let reordered = indices.iter().map(|&i| values[i].clone()).collect();
            *values = reordered;
        }
    }

    /// Apply a reordered reaction view back to the struct when the caller wants to persist sorting.
    fn apply_reordered_reaction_views(
        &mut self,
        views: ReorderedReactionViews,
        persist: bool,
    ) -> KineticsResult<()> {
        if !persist {
            return Ok(());
        }

        self.assign_reordered_reaction_views(views);
        self.finalize_reaction_mutation(true, None)?;
        Ok(())
    }

    /// Derive the best lifecycle state from the currently available data.
    fn state_from_available_data(&self) -> KinDataState {
        self.workflow_snapshot().visible_state()
    }

    /// Convert the structural branch classification into a workflow state.
    fn state_from_branch_kind(branch_kind: ReactionBranchKind) -> KinDataState {
        match branch_kind {
            ReactionBranchKind::Mechanism => KinDataState::MechanismBuilt,
            ReactionBranchKind::Parsed => KinDataState::ReactionDataParsed,
            ReactionBranchKind::Direct => KinDataState::DirectReactionsLoaded,
            ReactionBranchKind::Library => KinDataState::LibraryResolved,
            ReactionBranchKind::Shortcut => KinDataState::ShortcutsSelected,
            ReactionBranchKind::Empty => KinDataState::Empty,
            ReactionBranchKind::Normalized => KinDataState::Sorted,
        }
    }

    /// Check whether the computed state is actually compatible with the current payload.
    pub(crate) fn validate_state_contract(&self, state: KinDataState) -> KineticsResult<()> {
        let mismatch = |message: String| KineticsError::InvalidState(message);
        let snapshot = self.workflow_snapshot();

        match state {
            KinDataState::Empty => Ok(()),
            KinDataState::ShortcutsSelected
                if snapshot.branch_kind == ReactionBranchKind::Shortcut =>
            {
                Ok(())
            }
            KinDataState::LibraryResolved
                if snapshot.branch_kind == ReactionBranchKind::Library =>
            {
                Ok(())
            }
            KinDataState::DirectReactionsLoaded
                if snapshot.branch_kind == ReactionBranchKind::Direct =>
            {
                Ok(())
            }
            KinDataState::ReactionDataParsed if snapshot.has_parsed_branch => Ok(()),
            KinDataState::MechanismBuilt if snapshot.has_mechanism_branch => Ok(()),
            KinDataState::Analyzed
                if snapshot.has_analyzed_artifacts
                    || snapshot.has_parsed_branch
                    || snapshot.has_equations =>
            {
                Ok(())
            }
            KinDataState::Sorted if snapshot.normalized_view_ready => Ok(()),
            KinDataState::ShortcutsSelected => Err(mismatch(
                "ShortcutsSelected requires shortcut metadata or shortcut reaction ids".to_string(),
            )),
            KinDataState::LibraryResolved => Err(mismatch(
                "LibraryResolved requires library pairs or library reaction ids".to_string(),
            )),
            KinDataState::DirectReactionsLoaded => Err(mismatch(
                "DirectReactionsLoaded requires raw reaction values or document ids".to_string(),
            )),
            KinDataState::ReactionDataParsed => Err(mismatch(
                "ReactionDataParsed requires parsed reaction data".to_string(),
            )),
            KinDataState::MechanismBuilt => Err(mismatch(
                "MechanismBuilt requires reaction data from mechanism construction".to_string(),
            )),
            KinDataState::Analyzed => {
                if snapshot.has_analyzed_artifacts
                    || snapshot.has_parsed_branch
                    || snapshot.has_equations
                {
                    Ok(())
                } else {
                    Err(mismatch(
                        "Analyzed requires parsed reaction data, equations, or analyzed artifacts"
                            .to_string(),
                    ))
                }
            }
            KinDataState::Sorted => Err(mismatch(
                "Sorted requires a coherent normalized reaction view".to_string(),
            )),
        }
    }

    /// Recompute the visible lifecycle state after a structural mutation.
    pub(crate) fn refresh_state_from_available_data(&mut self) -> KineticsResult<()> {
        self.rebuild_reaction_ids_from_available_metadata();
        self.publish_refreshed_state(None)
    }

    /// Rebuild canonical ids from richer metadata when the current ids are missing or stale.
    fn rebuild_reaction_ids_from_available_metadata(&mut self) {
        let Some(count) = self.available_reaction_slot_count() else {
            return;
        };

        if self
            .canonical_reaction_ids()
            .is_some_and(|ids| ids.len() == count)
        {
            return;
        }

        let mut ids = Vec::with_capacity(count);
        for i in 0..count {
            ids.push(self.reaction_identity_from_metadata(i));
        }
        self.reaction_ids = Some(ids);
    }

    /// Report the widest reaction row count available across the current payload.
    fn available_reaction_slot_count(&self) -> Option<usize> {
        let count = [
            self.canonical_reaction_ids().map(|ids| ids.len()),
            self.reaction_data().map(|values| values.len()),
            self.reaction_values().map(|values| values.len()),
            self.shortcut_reactions.as_ref().map(|values| values.len()),
            self.vec_of_pairs.as_ref().map(|values| values.len()),
            Some(self.equations().len()),
            self.symbolic_constants().map(|values| values.len()),
            self.normalized_reactions().map(|values| values.len()),
        ]
        .into_iter()
        .flatten()
        .max();

        count.filter(|count| *count > 0)
    }

    /// Validate the canonical per-reaction contract before rebuild or export.
    fn validate_optional_len<T>(
        slot: Option<&[T]>,
        expected: usize,
        context: &'static str,
    ) -> KineticsResult<()> {
        if let Some(values) = slot {
            if values.len() != expected {
                return Err(KineticsError::LengthMismatch {
                    context,
                    left: values.len(),
                    right: expected,
                });
            }
        }
        Ok(())
    }

    fn validate_reaction_metadata_lengths(&self) -> KineticsResult<()> {
        let Some(vec_of_reaction_data) = self.reaction_data() else {
            if let Some(ids) = self.canonical_reaction_ids() {
                let expected = ids.len();
                Self::validate_optional_len(
                    Some(self.equations()),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(vec_of_equations vs reaction_ids)",
                )?;
                Self::validate_optional_len(
                    self.reaction_values(),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(vec_of_reaction_Values vs reaction_ids)",
                )?;
                Self::validate_optional_len(
                    self.shortcut_reactions.as_deref(),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(shortcut_reactions vs reaction_ids)",
                )?;
                Self::validate_optional_len(
                    self.vec_of_pairs.as_deref(),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(vec_of_pairs vs reaction_ids)",
                )?;
                Self::validate_optional_len(
                    self.symbolic_constants(),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(K_sym_vec vs reaction_ids)",
                )?;
                Self::validate_optional_len(
                    self.normalized_reactions(),
                    expected,
                    "KinData::validate_reaction_metadata_lengths(every_reaction vs reaction_ids)",
                )?;
            }
            return Ok(());
        };

        let expected = vec_of_reaction_data.len();
        Self::validate_optional_len(
            self.canonical_reaction_ids(),
            expected,
            "KinData::validate_reaction_metadata_lengths(reaction_ids vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            self.reaction_values(),
            expected,
            "KinData::validate_reaction_metadata_lengths(vec_of_reaction_Values vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            Some(self.equations()),
            expected,
            "KinData::validate_reaction_metadata_lengths(vec_of_equations vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            self.shortcut_reactions.as_deref(),
            expected,
            "KinData::validate_reaction_metadata_lengths(shortcut_reactions vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            self.vec_of_pairs.as_deref(),
            expected,
            "KinData::validate_reaction_metadata_lengths(vec_of_pairs vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            self.symbolic_constants(),
            expected,
            "KinData::validate_reaction_metadata_lengths(K_sym_vec vs vec_of_reaction_data)",
        )?;
        Self::validate_optional_len(
            self.normalized_reactions(),
            expected,
            "KinData::validate_reaction_metadata_lengths(every_reaction vs vec_of_reaction_data)",
        )?;
        Ok(())
    }

    /// Keep ids, derived views, and export views aligned after a mutation.
    fn sync_reaction_views_after_mutation(
        &mut self,
        rebuild_every_reaction: bool,
    ) -> KineticsResult<()> {
        self.rebuild_reaction_ids_from_available_metadata();
        let snapshot = self.capture_normalized_reaction_snapshot()?;

        if rebuild_every_reaction {
            if let Some(snapshot) = snapshot {
                self.every_reaction = Some(self.build_every_reaction_view(&snapshot)?);
            }
        } else {
            self.every_reaction = None;
        }

        Ok(())
    }

    /// Shared tail for mutation paths that need both derived views and visible state refreshed.
    ///
    /// When `state_override` is `Some`, the caller already knows the workflow stage that should
    /// be published after the cache update and we validate that stage explicitly.
    fn finalize_reaction_mutation(
        &mut self,
        rebuild_every_reaction: bool,
        state_override: Option<KinDataState>,
    ) -> KineticsResult<()> {
        self.sync_reaction_views_after_mutation(rebuild_every_reaction)?;
        self.publish_refreshed_state(state_override)
    }

    /// Publish the current visible state from the normalized payload or an explicit override.
    fn publish_refreshed_state(
        &mut self,
        state_override: Option<KinDataState>,
    ) -> KineticsResult<()> {
        let state = state_override.unwrap_or_else(|| self.state_from_available_data());
        self.set_state_checked(state)
    }

    /// Capture the current workflow shape once so state selection and validation use the same
    /// view of the payload.
    fn workflow_snapshot(&self) -> WorkflowSnapshot {
        WorkflowSnapshot {
            branch_kind: if self.normalized_reaction_view_is_ready() {
                ReactionBranchKind::Normalized
            } else if self.has_mechanism_branch() {
                ReactionBranchKind::Mechanism
            } else if self.has_parsed_branch() {
                ReactionBranchKind::Parsed
            } else if self.has_shortcut_branch() {
                ReactionBranchKind::Shortcut
            } else if self.has_library_branch() {
                ReactionBranchKind::Library
            } else if self.has_direct_branch() {
                ReactionBranchKind::Direct
            } else {
                ReactionBranchKind::Empty
            },
            normalized_view_ready: self.normalized_reaction_view_is_ready(),
            has_parsed_branch: self.has_parsed_branch(),
            has_equations: !self.equations().is_empty(),
            has_mechanism_branch: self.has_mechanism_branch(),
            has_analyzed_artifacts: self.has_analyzed_artifacts(),
        }
    }

    /// Report the canonical reaction count when enough data exists to derive it.
    fn canonical_reaction_count(&self) -> Option<usize> {
        // Canonical ids should drive the count whenever they are already available.
        self.canonical_reaction_ids()
            .map(|ids| ids.len())
            .or_else(|| self.reaction_data().map(|values| values.len()))
            .or_else(|| self.reaction_values().map(|values| values.len()))
    }

    /// Normalize the current reaction workspace before export or pretty-printing.
    fn normalize_before_use(&mut self) -> KineticsResult<()> {
        // Export and pretty-print operations should prefer the canonical id layer and the
        // parsed reaction rows, even when older cached views are stale.
        if self.normalized_reaction_view_is_ready() {
            return Ok(());
        }
        self.rebuild_every_reaction()?;
        Ok(())
    }

    /// Capture one coherent normalized reaction snapshot for rebuild and export paths.
    fn normalized_reaction_snapshot(
        &mut self,
    ) -> KineticsResult<Option<NormalizedReactionSnapshot>> {
        self.rebuild_reaction_ids_from_available_metadata();
        self.capture_normalized_reaction_snapshot()
    }

    /// Capture a coherent normalized reaction snapshot without mutating canonical ids first.
    ///
    /// Callers that already refreshed canonical ids can use this helper to avoid rebuilding the
    /// same ids twice in one mutation flow.
    fn capture_normalized_reaction_snapshot(
        &self,
    ) -> KineticsResult<Option<NormalizedReactionSnapshot>> {
        let Some(reaction_data) = self.reaction_data().map(|values| values.to_vec()) else {
            return Ok(None);
        };

        self.validate_reaction_metadata_lengths()?;
        let reaction_ids = self.canonical_reaction_ids().ok_or_else(|| {
            KineticsError::InvalidReactionData(
                "normalized reaction ids are not available".to_string(),
            )
        })?;

        Ok(Some(NormalizedReactionSnapshot {
            reaction_data,
            reaction_ids: reaction_ids.to_vec(),
            equations: self.equations().to_vec(),
            reaction_values: self.reaction_values().map(|values| values.to_vec()),
            k_sym: self.symbolic_constants().map(|values| values.to_vec()),
        }))
    }

    /// Rebuild the normalized `every_reaction` cache from canonical ids and parsed data.
    fn rebuild_every_reaction(&mut self) -> KineticsResult<()> {
        self.sync_reaction_views_after_mutation(true)
    }

    /// Build the export-ready reaction records from one coherent normalized snapshot.
    fn build_every_reaction_view(
        &self,
        snapshot: &NormalizedReactionSnapshot,
    ) -> KineticsResult<Vec<EveryReaction>> {
        // Keep the export view aligned with the canonical ids and the matching reaction row.
        if snapshot.reaction_ids.len() != snapshot.reaction_data.len() {
            return Err(KineticsError::LengthMismatch {
                context: "KinData::build_every_reaction_view(reaction_ids vs vec_of_reaction_data)",
                left: snapshot.reaction_ids.len(),
                right: snapshot.reaction_data.len(),
            });
        }
        snapshot
            .reaction_data
            .iter()
            .enumerate()
            .map(|(i, reaction)| {
                let equation = reaction.eq.clone();
                let (shortcut, lib_and_id) = match &snapshot.reaction_ids[i] {
                    ReactionIdentity::Shortcut(name) => (Some(name.clone()), None),
                    ReactionIdentity::LibraryReaction {
                        library,
                        reaction_id,
                    } => (None, Some((library.clone(), reaction_id.clone()))),
                    ReactionIdentity::Document { .. } => (None, None),
                };
                Ok(EveryReaction {
                    reaction_id: snapshot.reaction_ids[i].clone(),
                    shortcut,
                    lib_and_id,
                    reaction: reaction.clone(),
                    equation,
                    K_sym: snapshot.k_sym.as_ref().and_then(|v| v.get(i).cloned()),
                })
            })
            .collect()
    }
    /////////////////////////////////SETTING REACTIONS///////////////////////////////////////////
    /// Store direct equation strings as the active reaction set.
    pub fn set_reactions_directly(
        &mut self,
        reactions: Vec<String>,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> KineticsResult<()> {
        // Direct reactions are already the user-facing equation list, so store them as the
        // live equation cache and keep shortcut-specific fields empty.
        self.install_reaction_branch(ReactionBranchSnapshot::direct_equations(reactions))?;
        self.groups = groups;
        Ok(())
    }

    /// Store already parsed ReactionData rows as the active reaction set.
    ///
    /// This is the correct entry point for callers such as reactor helpers that construct
    /// ReactionData directly and should not manually write KinData internals.
    pub fn set_reaction_data_directly(
        &mut self,
        reaction_data: Vec<ReactionData>,
        groups: Option<HashMap<String, HashMap<String, usize>>>,
    ) -> KineticsResult<()> {
        self.install_reaction_branch(ReactionBranchSnapshot::parsed_reaction_data(reaction_data))?;
        self.groups = groups;
        self.refresh_substances_from_equations()?;
        self.publish_refreshed_state(Some(KinDataState::ReactionDataParsed))?;
        Ok(())
    }

    /// Expand a shortcut range like `C1..20` into canonical shortcut names.
    fn expand_shortcut_range(shortcut_range: &str) -> KineticsResult<Vec<String>> {
        let parts: Vec<&str> = shortcut_range.split("..").collect();
        if parts.len() != 2 {
            return Err(KineticsError::InvalidShortcutRange(
                shortcut_range.to_string(),
            ));
        }

        let prefix = parts[0]
            .chars()
            .take_while(|c| c.is_alphabetic())
            .collect::<String>();

        let start: usize = parts[0][prefix.len()..]
            .parse()
            .map_err(|_| KineticsError::InvalidShortcutRange(shortcut_range.to_string()))?;

        let end: usize = if let Some(last_char) = parts[1].chars().last() {
            if last_char.is_numeric() {
                parts[1]
                    .chars()
                    .rev()
                    .take_while(|c| c.is_numeric())
                    .collect::<String>()
                    .chars()
                    .rev()
                    .collect::<String>()
                    .parse()
                    .map_err(|_| KineticsError::InvalidShortcutRange(shortcut_range.to_string()))?
            } else {
                return Err(KineticsError::InvalidShortcutRange(
                    shortcut_range.to_string(),
                ));
            }
        } else {
            return Err(KineticsError::InvalidShortcutRange(
                shortcut_range.to_string(),
            ));
        };

        Ok((start..=end).map(|i| format!("{}_{}", prefix, i)).collect())
    }

    /// Expand a shortcut range like `C1..20` into `C_1`, `C_2`, ..., `C_20`.
    pub fn set_reactions_from_shortcut_range(
        &mut self,
        shortcut_range: String,
    ) -> KineticsResult<Vec<String>> {
        let res = Self::expand_shortcut_range(&shortcut_range)?;
        info!("task includes reactions as follows: {:#?}", &res);
        // Switching to a new shortcut branch should not keep raw payloads or parsed views from
        // the previous workflow.
        self.install_reaction_branch(ReactionBranchSnapshot::shortcut_selection(res.clone()))?;
        Ok(res)
    }

    /// decifer vector of shortcuts to full reaction names and store them in map {'library':"id of reaction"} and Vec<("library, reaction_id")>
    /// then open kinetic libraries, get reaction info corresponding to the guven shortcuts and store it in structure
    pub fn get_reactions_from_shortcuts(&mut self) -> KineticsResult<()> {
        let shortcut_reactions = if let Some(shortcuts) = self.shortcut_names() {
            shortcuts
        } else if self.canonical_reaction_ids().is_some() {
            return Err(KineticsError::InvalidReactionData(
                "canonical reaction ids are not shortcut ids".to_string(),
            ));
        } else {
            self.shortcut_reactions.clone().ok_or_else(|| {
                KineticsError::InvalidReactionData("no shortcut reactions available".to_string())
            })?
        };

        // Canonical shortcut ids win when they are available, so a stale legacy cache cannot
        // steer library resolution toward the wrong reactions.
        let vec: Vec<&str> = shortcut_reactions.iter().map(|s| s.as_str()).collect();
        let vec_of_pairs = decipher_vector_of_shortcuts_to_pairs(vec);
        // Cache loaded libraries so we only parse each JSON library once per workflow.
        let mut vec_of_reaction_values = Vec::new();
        let mut loaded_libraries: HashMap<String, KineticData> = HashMap::new();
        for (lib, reaction_id) in vec_of_pairs.iter() {
            let kin_instance = loaded_libraries
                .entry(lib.clone())
                .or_insert_with(KineticData::new);
            if kin_instance.LibKineticData.is_empty() || kin_instance.lib_name != *lib {
                kin_instance.open_json_files(lib)?;
            }
            let reaction_data_value = kin_instance.search_reactdata_by_reaction_id(reaction_id)?;
            vec_of_reaction_values.push(reaction_data_value);
        }
        self.install_reaction_branch(ReactionBranchSnapshot::library_resolved(
            Some(shortcut_reactions.to_vec()),
            vec_of_pairs.clone(),
            Some(vec_of_reaction_values),
        ))?;
        Ok(())
    }
    ///construct reaction mechanism
    pub fn construct_mechanism(
        &mut self,
        task_substances: Vec<String>,
        task_library: String,
    ) -> KineticsResult<()> {
        let mut found_mech = Mechanism_search::new(task_substances, task_library.clone());
        found_mech.mechfinder_api()?;
        // Mechanism construction replaces the current branch, so clear any prior raw payload.
        let vec_of_equations = found_mech.vec_of_reactions;
        if found_mech.reactdata.is_empty() {
            return Err(KineticsError::InvalidReactionData(
                "mechanism construction produced no reaction data".to_string(),
            ));
        }
        let vec_of_reaction_data = found_mech.reactdata;
        let reactions = found_mech.mechanism;
        self.install_reaction_branch(ReactionBranchSnapshot::mechanism_built(
            task_library,
            reactions,
            vec_of_reaction_data,
            vec_of_equations,
            found_mech.reactants,
        ))?;
        Ok(())
    }
    /////////////////////////////////REACTIONS MANIPULATIONS//////////////////////////////////////
    /// add manually reaction data as serde Value
    pub fn append_reaction(&mut self, reactions: Vec<Value>) -> KineticsResult<()> {
        let source = self.append_reaction_document_source();
        let values = if source.can_extend() {
            let mut combined = self
                .vec_of_reaction_Values
                .as_ref()
                .cloned()
                .unwrap_or_default();
            combined.extend(reactions);
            combined
        } else {
            reactions
        };
        self.install_reaction_branch(ReactionBranchSnapshot::raw_values(values, source, 0))
    }

    /// Decide whether a raw append can stay in the same document branch or must restart it.
    fn append_reaction_document_source(&self) -> DocumentSourceKind {
        let can_extend_raw_branch = self.reaction_ids.as_ref().is_none_or(|ids| {
            ids.iter().all(|reaction_id| {
                matches!(
                    reaction_id.document_source_kind(),
                    Some(DocumentSourceKind::Append | DocumentSourceKind::Raw)
                )
            })
        });

        if can_extend_raw_branch {
            DocumentSourceKind::Append
        } else {
            // A non-document branch cannot be extended in place, so the next raw payload starts fresh.
            DocumentSourceKind::Raw
        }
    }

    /// add manually reaction data as HashMaps
    pub fn append_reaction_from_map(
        &mut self,
        vec_of_maps: Vec<HashMap<String, Vec<f64>>>,
    ) -> KineticsResult<()> {
        use serde_json::json;
        let reactions: Vec<Value> = vec_of_maps.iter().map(|map| json!(map)).collect();
        self.append_reaction(reactions)
    }
    /// add manually reaction data with there shortcut names
    pub fn append_reaction_with_shortcut(
        &mut self,
        reactions: Vec<Value>,
        shortcuts: Vec<String>,
    ) -> KineticsResult<()> {
        if reactions.len() != shortcuts.len() {
            return Err(KineticsError::LengthMismatch {
                context: "KinData::append_reaction_with_shortcut",
                left: reactions.len(),
                right: shortcuts.len(),
            });
        }
        self.install_reaction_branch(ReactionBranchSnapshot::shortcut_values(
            shortcuts, reactions,
        ))
    }
    ///remove manually reaction data by its index
    pub fn remove_by_index(&mut self, index: usize) -> KineticsResult<()> {
        if index >= self.equations().len() {
            return Err(KineticsError::InvalidReactionData(format!(
                "reaction index {} is out of bounds for {} reactions",
                index,
                self.equations().len()
            )));
        }

        self.vec_of_equations.remove(index);
        Self::remove_option_index(&mut self.shortcut_reactions, index);
        Self::remove_option_index(&mut self.reaction_ids, index);
        Self::remove_option_index(&mut self.vec_of_reaction_Values, index);
        Self::remove_option_index(&mut self.vec_of_reaction_data, index);
        Self::remove_option_index(&mut self.vec_of_pairs, index);
        Self::remove_option_index(&mut self.K_sym_vec, index);
        Self::remove_option_index(&mut self.every_reaction, index);
        self.finalize_reaction_mutation(self.has_reaction_data(), None)?;
        Ok(())
    }
    pub fn remove_reaction_by_eq(&mut self, reaction_name: &str) -> KineticsResult<()> {
        let index = self
            .equations()
            .iter()
            .position(|eq_i| eq_i == reaction_name);
        if let Some(index) = index {
            self.remove_by_index(index)?;
            Ok(())
        } else {
            Err(KineticsError::InvalidReactionData(format!(
                "reaction equation `{}` was not found",
                reaction_name
            )))
        }
    }

    /////////////////////////////////COMPUTING AND PARSING REACTION DATA///////////////////////////////////////////
    /// parse serde Values with kinetic info into structs and store it in KinData structure
    pub fn reactdata_parsing(&mut self) -> KineticsResult<()> {
        self.ensure_state_allowed(
            "reaction-data parsing",
            &[
                KinDataState::Empty,
                KinDataState::ShortcutsSelected,
                KinDataState::LibraryResolved,
                KinDataState::DirectReactionsLoaded,
                KinDataState::MechanismBuilt,
                KinDataState::ReactionDataParsed,
                KinDataState::Analyzed,
                KinDataState::Sorted,
            ],
        )?;
        if let Some(vec_of_reaction_values) = self.reaction_values() {
            let (vec_ReactionData, vec_of_equations) =
                parse_kinetic_data_vec(vec_of_reaction_values.to_vec())?;
            //  info!("map_of_reaction_data: {:?}", &vec_of_reaction_data);
            //   info!("vec_of_equations: {:?}", &vec_of_equations);
            self.store_parsed_reaction_views(vec_ReactionData, vec_of_equations);
            self.finalize_reaction_mutation(true, Some(KinDataState::ReactionDataParsed))?;
            Ok(())
        } else {
            Err(KineticsError::InvalidReactionData(
                "no reaction values available".to_string(),
            ))
        }
    }
    // find reaction equations in reactdata
    pub fn equations_from_reactdata(&mut self) -> KineticsResult<()> {
        // Rebuild the equation cache from scratch to avoid duplicate entries.
        if let Some(vec_of_reaction_data) = self.reaction_data() {
            let vec_of_equations = vec_of_reaction_data
                .iter()
                .map(|reaction| reaction.eq.clone())
                .collect();
            self.store_parsed_reaction_views(vec_of_reaction_data.to_vec(), vec_of_equations);
            self.finalize_reaction_mutation(true, Some(KinDataState::ReactionDataParsed))?;
            Ok(())
        } else {
            Err(KineticsError::InvalidReactionData(
                "no reaction data available".to_string(),
            ))
        }
    }
    /// generates  Stoicheometric  data structures: matrix of stoicheometric coefficients, matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the
    /// kinetic function, G_matrix. As a rule, the degrees of concentration in the kinetic function coincide with the stoicheometric coefficients of
    /// the substances in the reaction; however, for empirical reactions they may differ from the stoicheometric coefficients.
    pub fn analyze_reactions(&mut self) -> KineticsResult<()> {
        self.ensure_state_allowed(
            "reaction analysis",
            &[
                KinDataState::DirectReactionsLoaded,
                KinDataState::ShortcutsSelected,
                KinDataState::LibraryResolved,
                KinDataState::MechanismBuilt,
                KinDataState::ReactionDataParsed,
                KinDataState::Analyzed,
                KinDataState::Sorted,
            ],
        )?;
        if self.equations().is_empty() {
            return Err(KineticsError::InvalidReactionData(
                "no reaction equations available".to_string(),
            ));
        }
        // iniciate instance of ReactionAnalyzer
        let mut StoichAnalyzer_instance = StoichAnalyzer::new();
        // copy vector of reactions to ReactionAnalyzer_instance
        StoichAnalyzer_instance.reactions = self.equations().to_vec();
        StoichAnalyzer_instance.groups = self.groups.clone();
        // Parse to find substance names. Commit the result only after the whole analysis
        // succeeds, so invalid formulae cannot leave a half-updated KinData workspace.
        let next_substances = if self.substances.is_empty() {
            // no need to parse for substances if they are already known
            StoichAnalyzer_instance.search_substances()?;
            StoichAnalyzer_instance.substances.clone()
        } else {
            self.substances.clone()
        };
        StoichAnalyzer_instance.substances = next_substances.clone();
        //find stoichiometric matrix and other matrices
        StoichAnalyzer_instance.analyse_reactions()?;
        StoichAnalyzer_instance.create_matrix_of_elements()?;
        self.substances = next_substances;
        self.stecheodata = StoichAnalyzer_instance;
        self.finalize_reaction_mutation(true, Some(KinDataState::Analyzed))?;
        Ok(())
    }
    /// parsing reaction data into structures and stoichometric calculations under one hood
    pub fn kinetic_main(&mut self) -> KineticsResult<()> {
        // Direct-reaction workflows already provide equations, so they can skip serde parsing.
        if self.reaction_values().is_some() {
            self.reactdata_parsing()?;
        } else if self.equations().is_empty() {
            return Err(KineticsError::InvalidReactionData(
                "no reaction values or equations available".to_string(),
            ));
        }
        self.analyze_reactions()?;
        Ok(())
    }
    ///////////////////////////INPUT/OUTPUT/////////////////////////////////////////////////////////
    /// save the chosen reaction kinetic data to json
    pub fn save_raw_reactions(&self, file_name: &str) -> KineticsResult<()> {
        let Some(vec_of_reaction_values) = self.reaction_values() else {
            return Err(KineticsError::InvalidReactionData(
                "no reaction values available".to_string(),
            ));
        };

        let json_array = json!(vec_of_reaction_values);
        let name = format!("{file_name}.json");
        let mut file = File::create(name)?;
        file.write_all(serde_json::to_string_pretty(&json_array)?.as_bytes())?;
        info!("Raw reactions have been written to raw_reactions.json");
        Ok(())
    }
    /// save the chosen reactions kinetic data with their shortcuts to json as dictionaries {"shortcut":"kinetic data"}
    pub fn save_reactions_with_shortcuts(&self, file_name: &str) -> KineticsResult<()> {
        let reaction_values = self.reaction_values().ok_or_else(|| {
            KineticsError::InvalidReactionData("no reaction values available".to_string())
        })?;
        let shortcuts = self.shortcut_names().ok_or_else(|| {
            KineticsError::InvalidReactionData("no vector of shortcuts available".to_string())
        })?;

        if shortcuts.len() != reaction_values.len() {
            return Err(KineticsError::LengthMismatch {
                context: "KinData::save_reactions_with_shortcuts",
                left: shortcuts.len(),
                right: reaction_values.len(),
            });
        }

        let map: HashMap<String, Value> = shortcuts
            .iter()
            .cloned()
            .zip(reaction_values.iter().cloned())
            .collect();

        let file = File::create(file_name)?;
        to_writer_pretty(file, &map)?;
        Ok(())
    }
    ///  algorithm: 1) takes the name of the file being created as an argument.
    ///  2) checks if there is a "KINETICS" or "REACTIONS" header there; if not, it adds them, and under this header, it adds records to a
    /// form available for use by serde (i.e., of the json type): "library_name":{"reaction_name":{ ...reaction data

    pub fn create_kinetics_document(
        &mut self,
        file_name: &str,
    ) -> KineticsResult<HashMap<String, HashMap<String, Value>>> {
        let path = Path::new(file_name);
        let file_exists = path.exists();

        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .append(true)
            .open(file_name)?;

        if !file_exists {
            writeln!(file, "KINETICS")?;
        } else {
            let reader = BufReader::new(File::open(file_name)?);
            let mut has_header = false;
            for line in reader.lines() {
                let line = line?;
                let normalized = line.trim().to_uppercase();
                if normalized == "KINETICS" || normalized == "REACTIONS" {
                    has_header = true;
                    break;
                }
            }

            if !has_header {
                writeln!(file, "\nKINETICS")?;
            }
        }
        let hashmap_to_save = self.kinetics_document_map()?;

        writeln!(file, "{}", serde_json::to_string_pretty(&hashmap_to_save)?)?;

        Ok(hashmap_to_save)
    }
    /// load kinetic data of reactions with their shortcuts from json, where they are stored as dictionaries {"shortcut":"kinetic data"}
    /// and save them into vector of kinetic data serde values and vector of shortcuts
    pub fn load_reactions_from_json(&mut self, file_name: &str) -> KineticsResult<()> {
        let file = File::open(file_name)?;
        let map: HashMap<String, Value> = from_reader(file)?;
        // Keep the load order deterministic so downstream exports and tests stay stable.
        let mut entries: Vec<(String, Value)> = map.into_iter().collect();
        entries.sort_by(|a, b| a.0.cmp(&b.0));
        let (shortcuts, reaction_values): (Vec<String>, Vec<Value>) = entries.into_iter().unzip();
        self.install_reaction_branch(ReactionBranchSnapshot::shortcut_values(
            shortcuts,
            reaction_values,
        ))?;

        Ok(())
    }

    /// Build the stoichiometric and elemental tables without printing them.
    pub fn substances_verbose_tables(&self) -> KineticsResult<(Table, Table)> {
        let subs = &self.substances;
        let reacts = self.equations();
        let sd = self.stoichiometric_analyzer();
        let elem_matrix = sd.matrix_of_elements.as_ref().ok_or_else(|| {
            KineticsError::InvalidReactionData(
                "stoichiometric element matrix is not available".to_string(),
            )
        })?;
        let unique_vec_of_elems = sd.unique_vec_of_elems.as_ref().ok_or_else(|| {
            KineticsError::InvalidReactionData("unique element list is not available".to_string())
        })?;
        let sm = sd.stecheo_matrx.clone();
        let n_react_from_matrix = sm.len();
        if n_react_from_matrix != reacts.len() {
            return Err(KineticsError::LengthMismatch {
                context: "stoichiometric matrix row count",
                left: n_react_from_matrix,
                right: reacts.len(),
            });
        }
        let n_subs_from_matrix = sm.first().map_or(0, |row| row.len());
        if n_subs_from_matrix != subs.len() {
            return Err(KineticsError::LengthMismatch {
                context: "stoichiometric matrix column count",
                left: n_subs_from_matrix,
                right: subs.len(),
            });
        }
        let (nrows, ncols) = elem_matrix.shape();
        if nrows != subs.len() {
            return Err(KineticsError::LengthMismatch {
                context: "element matrix row count",
                left: nrows,
                right: subs.len(),
            });
        }
        if ncols != unique_vec_of_elems.len() {
            return Err(KineticsError::LengthMismatch {
                context: "element matrix column count",
                left: ncols,
                right: unique_vec_of_elems.len(),
            });
        }

        let mut table_stecheo = Table::new();
        let mut header_row = vec![Cell::new("Reactions/Substances")];
        for sub in subs {
            header_row.push(Cell::new(sub));
        }
        table_stecheo.add_row(Row::new(header_row));

        for (i, react) in reacts.iter().enumerate() {
            let mut row = vec![Cell::new(react)];
            for j in 0..subs.len() {
                row.push(Cell::new(&format!("{:.4}", sm[i][j])));
            }
            table_stecheo.add_row(Row::new(row));
        }

        let mut elem_table = Table::new();
        let mut header_row = vec![Cell::new("Substances/Elements")];
        for elems in unique_vec_of_elems.iter() {
            header_row.push(Cell::new(elems));
        }
        elem_table.add_row(Row::new(header_row));
        for (i, sub) in subs.iter().enumerate() {
            let mut row = vec![Cell::new(sub)];
            for j in 0..unique_vec_of_elems.len() {
                let value = elem_matrix.get((i, j)).ok_or_else(|| {
                    KineticsError::InvalidReactionData(format!(
                        "missing element matrix value at row {}, column {}",
                        i, j
                    ))
                })?;
                row.push(Cell::new(&format!("{:.4}", value)));
            }
            elem_table.add_row(Row::new(row));
        }

        Ok((table_stecheo, elem_table))
    }

    /// Print the stoichiometric and elemental tables in a human-readable format.
    pub fn pretty_print_substances_verbose(&self) -> KineticsResult<()> {
        let (table_stecheo, elem_table) = self.substances_verbose_tables()?;
        table_stecheo.printstd();
        elem_table.printstd();
        Ok(())
    }

    /// Build the normalized reaction table without mutating the workflow state.
    pub fn normalized_reaction_data_table(&self) -> KineticsResult<Table> {
        self.build_reaction_data_table()
    }

    /// Normalize the workflow state if needed and then build the reaction table.
    pub fn reaction_data_table(&mut self) -> KineticsResult<Table> {
        self.normalize_before_use()?;
        self.normalized_reaction_data_table()
    }

    /// Render the reaction table without changing the workflow state.
    pub fn pretty_reaction_data_table(&self) -> KineticsResult<Table> {
        self.normalized_reaction_data_table()
    }

    /// Print the normalized reaction table in a human-readable format.
    pub fn pretty_print_reaction_data(&self) -> KineticsResult<()> {
        // Reuse the already-normalized cache when it is available; otherwise normalize a clone
        // so the live workspace stays untouched.
        let table = if self.normalized_reaction_view_is_ready() {
            self.pretty_reaction_data_table()?
        } else {
            let mut normalized = self.clone();
            normalized.normalize_before_use()?;
            normalized.pretty_reaction_data_table()?
        };
        table.printstd();
        Ok(())
    }

    /// Print the normalized reaction table without touching the workflow state.
    pub fn pretty_print_reaction_data_from_normalized(&self) -> KineticsResult<()> {
        let table = self.pretty_reaction_data_table()?;
        table.printstd();
        Ok(())
    }

    /// Build the normalized reaction table from the already normalized cache.
    fn build_reaction_data_table(&self) -> KineticsResult<Table> {
        let Some(every_reaction) = self.normalized_reactions() else {
            return Err(KineticsError::InvalidReactionData(
                "no normalized reaction view available".to_string(),
            ));
        };

        Self::build_reaction_data_table_from_entries(every_reaction)
    }

    /// Build a formatted reaction table from an already prepared normalized view.
    fn build_reaction_data_table_from_entries(
        every_reaction: &[EveryReaction],
    ) -> KineticsResult<Table> {
        let mut table = Table::new();
        table.add_row(row![
            "Reaction number",
            "Canonical ID",
            "Reaction Type",
            "Equation",
            "Reactants",
            "Kinetics Data"
        ]);

        for (i, reaction_entry) in every_reaction.iter().enumerate() {
            let reaction_type = format!("{:?}", reaction_entry.reaction.reaction_type);
            let equation = &reaction_entry.equation;
            let reactants = match &reaction_entry.reaction.react {
                Some(react) => format!("{:?}", react),
                None => "None".to_string(),
            };
            let kinetics_data = match &reaction_entry.reaction.data {
                ReactionKinetics::Elementary(data) => format!("{:?}", data),
                ReactionKinetics::Falloff(data) => format!("{:?}", data),
                ReactionKinetics::Pressure(data) => format!("{:?}", data),
                ReactionKinetics::ThreeBody(data) => format!("{:?}", data),
            };

            table.add_row(Row::new(vec![
                Cell::new(i.to_string().as_str()),
                Cell::new(&format!("{:?}", reaction_entry.reaction_id)),
                Cell::new(&reaction_type),
                Cell::new(equation),
                Cell::new(&reactants),
                Cell::new(&kinetics_data),
            ]));
        }

        Ok(table)
    }

    /// Build the kinetics document payload without writing it to disk.
    pub fn kinetics_document_map(
        &mut self,
    ) -> KineticsResult<HashMap<String, HashMap<String, Value>>> {
        self.normalize_before_use()?;
        self.kinetics_document_map_from_normalized()
    }

    /// Build the kinetics export payload from the already normalized cache.
    pub fn kinetics_document_map_from_normalized(
        &self,
    ) -> KineticsResult<HashMap<String, HashMap<String, Value>>> {
        let every_reaction = self.normalized_reactions().ok_or_else(|| {
            KineticsError::InvalidReactionData("no normalized reaction view available".to_string())
        })?;
        // The normalized export view is the contract here; stale raw payload caches must not
        // block export once the coherent per-reaction view already exists.
        if every_reaction.is_empty() {
            // Keep the empty-view error aligned with the broader normalized-view contract.
            return Err(KineticsError::InvalidReactionData(
                "no normalized reaction view available".to_string(),
            ));
        }

        let mut hashmap_to_save: HashMap<String, HashMap<String, Value>> = HashMap::new();
        for reaction in every_reaction {
            let (library_name, reaction_name) =
                Self::reaction_id_to_document_key(&reaction.reaction_id);
            let reactdata = reaction.reaction.to_serde_value()?;
            hashmap_to_save
                .entry(library_name)
                .or_default()
                .insert(reaction_name, reactdata);
        }

        Ok(hashmap_to_save)
    }
    /// both reaction data and matrces are pretty printed
    pub fn pretty_print_kindata(&self) -> KineticsResult<()> {
        self.pretty_print_reaction_data()?;
        self.pretty_print_substances_verbose()?;
        Ok(())
    }

    /// Print the combined output without mutating the workflow state.
    pub fn pretty_print_kindata_from_normalized(&self) -> KineticsResult<()> {
        self.pretty_print_reaction_data_from_normalized()?;
        self.pretty_print_substances_verbose()?;
        Ok(())
    }

    /////////////////////KINETIC CONTANT FUNCTIONS //////////////////////
    ///  Returns the kinetic constant for a reaction with number reaction_id
    pub fn calc_K_const_for_1_react(
        &self,
        reaction_id: usize,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> KineticsResult<f64> {
        let Some(reactions) = self.reaction_data() else {
            return Err(KineticsError::InvalidReactionData(
                "no vector of reaction data available".to_string(),
            ));
        };

        let reaction_data = reactions.get(reaction_id).ok_or_else(|| {
            KineticsError::InvalidReactionData(format!(
                "reaction index {} is out of bounds for {} reactions",
                reaction_id,
                reactions.len()
            ))
        })?;
        reaction_data.K_const(temp, pres, concentrations)
    }
    /// 1) Returns the kinetic constant for all reactions found 2) sort reactions by reaction rate constant and regroup  3)
    /// if  save_rearranged = Some(true) rearrange all data in structure according to reaction rate constant magnitude
    /// 4) pretty print table with headera    "Reaction number",  "Equation",  "log10 of Constant at given T"
    pub fn calc_K_const_for_all_reactions(
        &mut self,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
        save_rearranged: Option<bool>,
    ) -> KineticsResult<Vec<f64>> {
        let Some(snapshot) = self.normalized_reaction_snapshot()? else {
            return Err(KineticsError::InvalidReactionData(
                "no vector of reaction data available".to_string(),
            ));
        };
        let reactions = &snapshot.reaction_data;

        let mut vec_of_K_const = Vec::new();
        for reaction in reactions.iter() {
            let K = reaction.K_const(temp, pres, concentrations.clone())?;
            vec_of_K_const.push(K);
        }

        let mut indices: Vec<usize> = (0..vec_of_K_const.len()).collect();
        indices.sort_by(|&i, &j| {
            vec_of_K_const[j]
                .partial_cmp(&vec_of_K_const[i])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        vec_of_K_const = indices.iter().map(|&i| vec_of_K_const[i]).collect();
        let reordered_views = snapshot.reordered(&indices)?;
        self.apply_reordered_reaction_views(reordered_views, save_rearranged == Some(true))?;
        Ok(vec_of_K_const)
    }
    /// the same but for every reaction for given range of T calculated K values and maximum values compared to those for another reaction
    pub fn calc_K_const_for_all_reactions_forTrange(
        &mut self,
        T0: f64,
        Tend: f64,
        n: usize,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
        save_rearranged: Option<bool>,
    ) -> KineticsResult<Vec<f64>> {
        let Some(snapshot) = self.normalized_reaction_snapshot()? else {
            return Err(KineticsError::InvalidReactionData(
                "no vector of reaction data available".to_string(),
            ));
        };
        let reactions = &snapshot.reaction_data;

        let mut vec_of_K_const = Vec::new();
        for reaction in reactions.iter() {
            let K_vec = reaction.K_const_for_T_range(T0, Tend, n, pres, concentrations.clone())?;

            let Some(max_K) = K_vec.iter().max_by(|a, b| a.total_cmp(b)) else {
                return Err(KineticsError::InvalidReactionData(
                    "temperature range produced no rate constants".to_string(),
                ));
            };
            vec_of_K_const.push(*max_K);
        }

        let mut indices: Vec<usize> = (0..vec_of_K_const.len()).collect();
        indices.sort_by(|&i, &j| {
            vec_of_K_const[j]
                .partial_cmp(&vec_of_K_const[i])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        vec_of_K_const = indices.iter().map(|&i| vec_of_K_const[i]).collect();
        let reordered_views = snapshot.reordered(&indices)?;
        self.apply_reordered_reaction_views(reordered_views, save_rearranged == Some(true))?;
        Ok(vec_of_K_const)
    }
    pub fn calc_sym_constants(
        &mut self,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, Expr>>,
        T_scaling: Option<Expr>,
    ) -> KineticsResult<()> {
        self.ensure_state_allowed(
            "symbolic constant calculation",
            &[
                KinDataState::Empty,
                KinDataState::Analyzed,
                KinDataState::Sorted,
            ],
        )?;
        let Some(vec_of_reaction_data) = self.reaction_data() else {
            return Err(KineticsError::InvalidReactionData(
                "no reaction data available".to_string(),
            ));
        };
        let mut vec_of_k_sym: Vec<Expr> = Vec::new();
        for reaction in vec_of_reaction_data.iter() {
            let K_sym =
                reaction.K_sym_with_scaled_T(pres, concentrations.clone(), T_scaling.clone())?;
            vec_of_k_sym.push(K_sym);
        }
        self.store_symbolic_constants(vec_of_k_sym);
        self.finalize_reaction_mutation(true, Some(KinDataState::Analyzed))?;
        Ok(())
    }
}
