//! # Molecular Mass and Atomic Composition Module
//!
//! ## Purpose
//! This module provides comprehensive functionality for parsing chemical formulas and calculating
//! molecular properties. It handles complex chemical notation including parentheses, stoichiometric
//! coefficients, phase markers, and custom chemical groups.
//!
//! ## Main Data Structures
//! - `Element`: Structure containing element name and atomic mass
//! - `ELEMENTS`: Static array of periodic table elements with atomic masses
//! - Chemical formula parser that returns `HashMap<String, usize>` (element -> count)
//! - Matrix operations using nalgebra::DMatrix for elemental composition matrices
//!
//! ## Key Logic Implementation
//! 1. **Formula Parsing**: Complex regex-based parser handling nested parentheses and coefficients
//! 2. **Phase Filtering**: Removes phase markers like (g), (l), (s), (c) from formulas
//! 3. **Group Expansion**: Converts chemical groups (e.g., "Me" -> {"C":1, "H":3}) into elements
//! 4. **Bracket Processing**: Handles nested parentheses with multipliers like Ca(NO3)2
//! 5. **Matrix Construction**: Creates elemental composition matrices for multiple substances
//!
//! ## Usage Patterns
//! ```rust, ignore
//! // Basic usage
//! let (mass, composition) = calculate_molar_mass("C6H8O6".to_string(), None).unwrap();
//!
//! // With chemical groups
//! let groups = Some(HashMap::from([("Me".to_string(),
//!     HashMap::from([("C".to_string(), 1), ("H".to_string(), 3)]))]);
//! let composition = parse_formula("C6H5Me".to_string(), groups).unwrap();
//!
//! // Matrix operations
//! let (matrix, elements) = create_elem_composition_matrix(vec!["H2O", "CO2"], None).unwrap();
//! ```
//!
//! ## Interesting Features
//! - **Intelligent Parsing**: Handles complex cases like "Ca(NO3)2", "C(OOH)2N(ClO)3"
//! - **Phase Awareness**: Automatically strips phase indicators from formulas
//! - **Group Support**: Extensible system for custom chemical groups (methyl, phenyl, etc.)
//! - **Matrix Generation**: Creates elemental composition matrices for stoichiometric analysis
//! - **Robust Error Handling**: Comprehensive validation and error reporting
//! - **Case Sensitivity**: Smart handling of element capitalization (Ca vs CA)

use crate::Kinetics::error::{KineticsError, KineticsResult};
use nalgebra::DMatrix;
use std::collections::{BTreeSet, HashMap};

const MAX_FORMULA_BYTES: usize = 4096;
const MAX_FORMULA_NESTING: usize = 64;

// Define a struct to hold element data
pub struct Element {
    name: &'static str,
    atomic_mass: f64,
}

// Define a list of elements and their atomic masses
pub const ELEMENTS: &[Element] = &[
    Element {
        name: "H",
        atomic_mass: 1.008,
    },
    Element {
        name: "He",
        atomic_mass: 4.0026,
    },
    Element {
        name: "Li",
        atomic_mass: 6.94,
    },
    Element {
        name: "Be",
        atomic_mass: 9.0122,
    },
    Element {
        name: "B",
        atomic_mass: 10.81,
    },
    Element {
        name: "C",
        atomic_mass: 12.011,
    },
    Element {
        name: "N",
        atomic_mass: 14.007,
    },
    Element {
        name: "O",
        atomic_mass: 15.999,
    },
    Element {
        name: "F",
        atomic_mass: 18.998,
    },
    Element {
        name: "Ne",
        atomic_mass: 20.18,
    },
    Element {
        name: "Na",
        atomic_mass: 22.99,
    },
    Element {
        name: "Mg",
        atomic_mass: 24.305,
    },
    Element {
        name: "Al",
        atomic_mass: 26.98,
    },
    Element {
        name: "Si",
        atomic_mass: 28.085,
    },
    Element {
        name: "P",
        atomic_mass: 30.974,
    },
    Element {
        name: "S",
        atomic_mass: 32.065,
    },
    Element {
        name: "Cl",
        atomic_mass: 35.45,
    },
    Element {
        name: "Ar",
        atomic_mass: 39.948,
    },
    Element {
        name: "K",
        atomic_mass: 39.102,
    },
    Element {
        name: "Ca",
        atomic_mass: 40.08,
    },
    Element {
        name: "Sc",
        atomic_mass: 44.9559,
    },
    Element {
        name: "Ti",
        atomic_mass: 47.867,
    },
    Element {
        name: "V",
        atomic_mass: 50.9415,
    },
    Element {
        name: "Cr",
        atomic_mass: 51.9961,
    },
    Element {
        name: "Mn",
        atomic_mass: 54.938,
    },
    Element {
        name: "Fe",
        atomic_mass: 55.845,
    },
    Element {
        name: "Co",
        atomic_mass: 58.933,
    },
    Element {
        name: "Ni",
        atomic_mass: 58.69,
    },
    Element {
        name: "Cu",
        atomic_mass: 63.546,
    },
    Element {
        name: "Zn",
        atomic_mass: 65.38,
    },
    Element {
        name: "Ga",
        atomic_mass: 69.723,
    },
    Element {
        name: "Ge",
        atomic_mass: 72.64,
    },
    Element {
        name: "As",
        atomic_mass: 74.9216,
    },
    Element {
        name: "Se",
        atomic_mass: 78.96,
    },
    Element {
        name: "Br",
        atomic_mass: 79.904,
    },
    Element {
        name: "Kr",
        atomic_mass: 83.798,
    },
    Element {
        name: "Rb",
        atomic_mass: 85.4678,
    },
    Element {
        name: "Sr",
        atomic_mass: 87.62,
    },
    Element {
        name: "Y",
        atomic_mass: 88.9059,
    },
    Element {
        name: "Zr",
        atomic_mass: 91.224,
    },
    Element {
        name: "Nb",
        atomic_mass: 92.9064,
    },
    Element {
        name: "Mo",
        atomic_mass: 95.94,
    },
    Element {
        name: "Tc",
        atomic_mass: 98.0,
    },
    Element {
        name: "Ru",
        atomic_mass: 101.07,
    },
    // Add more elements here...
    // Add more elements here...
];

/// Returns whether `symbol` is a known chemical element in the parser's
/// periodic-table data.
///
/// Formula parsing intentionally accepts syntactically valid uppercase tokens
/// before an application decides what they mean. Callers which construct a
/// physical element inventory can use this predicate to reject an unknown
/// token without maintaining a second periodic-table list.
pub fn is_known_element(symbol: &str) -> bool {
    ELEMENTS.iter().any(|element| element.name == symbol)
}

fn strip_phase_suffix(formula: &str) -> &str {
    ["(g)", "(G)", "(l)", "(L)", "(s)", "(S)", "(c)", "(C)"]
        .into_iter()
        .find_map(|suffix| formula.strip_suffix(suffix))
        .unwrap_or(formula)
}
// Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
// so we need to convert them into regular elements
fn handle_groups(
    counts: HashMap<String, usize>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
    formula: &str,
) -> KineticsResult<HashMap<String, usize>> {
    let Some(groups) = groups else {
        return Ok(counts);
    };
    let mut expanded = HashMap::new();
    for (symbol, count) in counts {
        let Some(composition) = groups.get(&symbol) else {
            let total = expanded.entry(symbol.clone()).or_insert(0usize);
            *total = total.checked_add(count).ok_or_else(|| {
                formula_error(formula, format!("atom count overflow for `{symbol}`"))
            })?;
            continue;
        };
        if symbol.is_empty()
            || !symbol.is_ascii()
            || !symbol.chars().all(char::is_alphabetic)
            || is_known_element(&symbol)
        {
            return Err(KineticsError::InvalidCustomGroup {
                group: symbol,
                reason: "group name must be alphabetic and must not shadow an element".into(),
            });
        }
        for (element, group_count) in composition {
            if groups.contains_key(element) {
                return Err(KineticsError::InvalidCustomGroup {
                    group: symbol,
                    reason: format!("nested group `{element}` is unsupported"),
                });
            }
            if !is_known_element(element) {
                return Err(KineticsError::UnknownElement {
                    formula: formula.to_string(),
                    symbol: element.clone(),
                });
            }
            let amount = group_count.checked_mul(count).ok_or_else(|| {
                formula_error(formula, format!("atom count overflow for `{element}`"))
            })?;
            let total = expanded.entry(element.clone()).or_insert(0usize);
            *total = total.checked_add(amount).ok_or_else(|| {
                formula_error(formula, format!("atom count overflow for `{element}`"))
            })?;
        }
    }
    Ok(expanded)
}

fn validate_custom_groups(
    groups: &HashMap<String, HashMap<String, usize>>,
    formula: &str,
) -> KineticsResult<()> {
    let mut names: Vec<_> = groups.keys().collect();
    names.sort();
    for name in names.iter() {
        if name.is_empty()
            || !name.is_ascii()
            || !name.chars().all(char::is_alphabetic)
            || is_known_element(name)
        {
            return Err(KineticsError::InvalidCustomGroup {
                group: (*name).clone(),
                reason: "group name must be ASCII alphabetic and must not shadow an element".into(),
            });
        }
        let composition = &groups[*name];
        if composition.is_empty() {
            return Err(KineticsError::InvalidCustomGroup {
                group: (*name).clone(),
                reason: "composition must contain at least one element".into(),
            });
        }
        for element in composition.keys() {
            if groups.contains_key(element) {
                return Err(KineticsError::InvalidCustomGroup {
                    group: (*name).clone(),
                    reason: format!("nested or cyclic group `{element}` is unsupported"),
                });
            }
            if !is_known_element(element) {
                return Err(KineticsError::UnknownElement {
                    formula: formula.to_string(),
                    symbol: element.clone(),
                });
            }
            if composition[element] == 0 {
                return Err(KineticsError::InvalidCustomGroup {
                    group: (*name).clone(),
                    reason: format!("atom count for `{element}` must be positive"),
                });
            }
        }
    }
    Ok(())
}

fn validate_exact_group_composition(
    group: &str,
    composition: &HashMap<String, usize>,
    groups: &HashMap<String, HashMap<String, usize>>,
    formula: &str,
) -> KineticsResult<()> {
    if composition.is_empty() {
        return Err(KineticsError::InvalidCustomGroup {
            group: group.to_string(),
            reason: "composition must contain at least one element".into(),
        });
    }
    for (element, count) in composition {
        if *count == 0 {
            return Err(KineticsError::InvalidCustomGroup {
                group: group.to_string(),
                reason: format!("atom count for `{element}` must be positive"),
            });
        }
        if !is_known_element(element) {
            if groups.contains_key(element) {
                return Err(KineticsError::InvalidCustomGroup {
                    group: group.to_string(),
                    reason: format!("nested or cyclic group `{element}` is unsupported"),
                });
            }
            return Err(KineticsError::UnknownElement {
                formula: formula.to_string(),
                symbol: element.clone(),
            });
        }
    }
    Ok(())
}
fn formula_error(formula: &str, message: impl Into<String>) -> KineticsError {
    KineticsError::FormulaSyntax {
        formula: formula.to_string(),
        reason: message.into(),
    }
}

fn char_at(formula: &str, idx: usize) -> KineticsResult<char> {
    formula
        .as_bytes()
        .get(idx)
        .map(|byte| *byte as char)
        .ok_or_else(|| formula_error(formula, format!("missing character at position {}", idx)))
}

fn parse_usize_slice(formula: &str, start: usize, end: usize) -> KineticsResult<usize> {
    formula[start..end].parse::<usize>().map_err(|_| {
        formula_error(
            formula,
            format!("invalid integer slice `{}`", &formula[start..end]),
        )
    })
}
/// Parse a numeric suffix that starts at `start` and return both the value and the first
/// index after the suffix.
fn parse_stoichio_suffix(formula: &str, start: usize) -> KineticsResult<(usize, usize)> {
    let mut end = start + 1;
    while end < formula.len() && char_at(formula, end)?.is_digit(10) {
        end += 1;
    }
    let count = parse_usize_slice(formula, start, end)?;
    if count == 0 {
        return Err(formula_error(
            formula,
            format!("stoichiometric count at position {start} must be positive"),
        ));
    }
    Ok((count, end))
}

fn merge_counts(
    target: &mut HashMap<String, usize>,
    source: HashMap<String, usize>,
    multiplier: usize,
    formula: &str,
) -> KineticsResult<()> {
    for (atom, count) in source {
        let amount = count.checked_mul(multiplier).ok_or_else(|| {
            formula_error(
                formula,
                format!("atom count multiplication overflow for `{atom}`"),
            )
        })?;
        let total = target.entry(atom.clone()).or_insert(0usize);
        *total = total.checked_add(amount).ok_or_else(|| {
            formula_error(
                formula,
                format!("atom count addition overflow for `{atom}`"),
            )
        })?;
    }
    Ok(())
}

/// Parse one formula segment until a closing parenthesis or the end of the string.
fn parse_formula_segment(
    formula: &str,
    mut i: usize,
    group_names: &[String],
    in_brackets: bool,
    initial_formula: &str,
    depth: usize,
) -> KineticsResult<(HashMap<String, usize>, usize)> {
    if depth > MAX_FORMULA_NESTING {
        return Err(KineticsError::FormulaLimitExceeded {
            resource: "nesting depth",
            requested: depth,
            limit: MAX_FORMULA_NESTING,
        });
    }
    let mut counts = HashMap::new();

    while i < formula.len() {
        let ch = char_at(formula, i)?;

        if ch == ')' {
            if in_brackets {
                return Ok((counts, i + 1));
            }
            return Err(formula_error(
                initial_formula,
                "unexpected closing parenthesis",
            ));
        }

        if ch == '(' {
            let (inner_counts, next_i) = parse_formula_segment(
                formula,
                i + 1,
                group_names,
                true,
                initial_formula,
                depth + 1,
            )?;
            i = next_i;

            let multiplier = if i < formula.len() && char_at(formula, i)?.is_digit(10) {
                let (multiplier, next_after_multiplier) = parse_stoichio_suffix(formula, i)?;
                i = next_after_multiplier;
                multiplier
            } else {
                1
            };

            merge_counts(&mut counts, inner_counts, multiplier, initial_formula)?;
            continue;
        }

        if let Some(group_name) = group_names
            .iter()
            .find(|group_name| formula[i..].starts_with(group_name.as_str()))
        {
            i += group_name.len();
            let multiplier = if i < formula.len() && char_at(formula, i)?.is_digit(10) {
                let (multiplier, next_after_multiplier) = parse_stoichio_suffix(formula, i)?;
                i = next_after_multiplier;
                multiplier
            } else {
                1
            };

            let total = counts.entry(group_name.clone()).or_insert(0usize);
            *total = total.checked_add(multiplier).ok_or_else(|| {
                formula_error(
                    initial_formula,
                    format!("count overflow for `{group_name}`"),
                )
            })?;
            continue;
        }

        if ch.is_uppercase() {
            let start = i;
            i += 1;
            if i < formula.len() {
                let next_char = char_at(formula, i)?;
                if next_char.is_lowercase() {
                    i += 1;
                }
            }

            let element = &formula[start..i];
            let multiplier = if i < formula.len() && char_at(formula, i)?.is_digit(10) {
                let (multiplier, next_after_multiplier) = parse_stoichio_suffix(formula, i)?;
                i = next_after_multiplier;
                multiplier
            } else {
                1
            };

            let total = counts.entry(element.to_string()).or_insert(0usize);
            *total = total.checked_add(multiplier).ok_or_else(|| {
                formula_error(initial_formula, format!("count overflow for `{element}`"))
            })?;
            continue;
        }

        if ch.is_digit(10) {
            return Err(formula_error(
                initial_formula,
                format!("unexpected stoichiometric coefficient at position {}", i),
            ));
        }

        return Err(formula_error(
            initial_formula,
            format!("failed to parse a token at position {}", i),
        ));
    }

    if in_brackets {
        return Err(formula_error(
            initial_formula,
            "missing closing parenthesis",
        ));
    }

    Ok((counts, i))
}

/// Parse a formula or an exact species key with an explicitly supplied composition. Exact keys
/// are resolved before tokenization; otherwise `groups` supplies abbreviation tokens such as
/// Me/Ph, matched longest-first when names share a prefix.
pub fn parse_formula(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<HashMap<String, usize>> {
    if !formula.is_ascii() {
        return Err(formula_error(&formula, "formula syntax must be ASCII"));
    }
    if formula.len() > MAX_FORMULA_BYTES {
        return Err(KineticsError::FormulaLimitExceeded {
            resource: "byte length",
            requested: formula.len(),
            limit: MAX_FORMULA_BYTES,
        });
    }
    let initial_formula = formula.clone();
    let formula = strip_phase_suffix(&formula.replace(' ', "")).to_string();
    if formula.is_empty() {
        return Err(formula_error(
            &initial_formula,
            "formula must contain at least one atom or group",
        ));
    }
    if let Some(groups) = groups.as_ref() {
        if let Some(composition) = groups.get(&formula) {
            validate_exact_group_composition(&formula, composition, groups, &initial_formula)?;
            return Ok(composition.clone());
        }
        validate_custom_groups(groups, &initial_formula)?;
    }

    let mut group_names = groups
        .as_ref()
        .map(|group_map| {
            let mut names = group_map.keys().cloned().collect::<Vec<_>>();
            names.sort_by(|a, b| b.len().cmp(&a.len()));
            names
        })
        .unwrap_or_default();

    // Element symbols remain available through uppercase/lowercase matching. Custom groups are
    // matched longest-first so prefix-related aliases have deterministic tokenization.
    if group_names.is_empty() {
        group_names = Vec::new();
    }

    let (counts, _) = parse_formula_segment(&formula, 0, &group_names, false, &initial_formula, 0)?;
    handle_groups(counts, groups, &initial_formula)
} // end of parse_formula

// Function to calculate the molar mass of a substance given its chemical formula
pub fn calculate_molar_mass(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<(f64, HashMap<String, usize>)> {
    let counts = parse_formula(formula.clone(), groups)?;
    let molar_mass = molar_mass_for_composition(counts.clone(), &formula)?;
    Ok((molar_mass, counts))
}

fn molar_mass_for_composition(
    counts: HashMap<String, usize>,
    formula: &str,
) -> KineticsResult<f64> {
    if counts.is_empty() {
        return Err(formula_error(
            formula,
            "composition must contain at least one atom",
        ));
    }
    let mut molar_mass = 0.0;
    for (element, count) in counts {
        if count == 0 {
            return Err(formula_error(
                formula,
                format!("atom count for `{element}` must be positive"),
            ));
        }
        let atomic_mass = ELEMENTS
            .iter()
            .find(|record| record.name == element)
            .ok_or_else(|| KineticsError::UnknownElement {
                formula: formula.to_string(),
                symbol: element.clone(),
            })?
            .atomic_mass;
        molar_mass += atomic_mass * count as f64;
    }
    if !molar_mass.is_finite() {
        return Err(formula_error(formula, "computed molar mass is non-finite"));
    }
    Ok(molar_mass)
}

pub fn calculate_molar_mass_for_composition(counts: HashMap<String, usize>) -> KineticsResult<f64> {
    molar_mass_for_composition(counts, "<element composition>")
}

fn validate_physical_composition(
    counts: &HashMap<String, usize>,
    formula: &str,
) -> KineticsResult<()> {
    if counts.is_empty() {
        return Err(formula_error(
            formula,
            "composition must contain at least one atom",
        ));
    }
    for symbol in counts.keys() {
        if counts[symbol] == 0 {
            return Err(formula_error(
                formula,
                format!("atom count for `{symbol}` must be positive"),
            ));
        }
        if !is_known_element(symbol) {
            return Err(KineticsError::UnknownElement {
                formula: formula.to_string(),
                symbol: symbol.clone(),
            });
        }
    }
    Ok(())
}

// Function to calculate the molar mass of a vector of chemical formulas
pub fn calculate_molar_mass_of_vector_of_subs(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<Vec<f64>> {
    let mut molar_masses = Vec::new();
    for formula in vec_of_formulae.iter() {
        let (molar_mass, _) = calculate_molar_mass((*formula).to_string(), groups.clone())?;
        molar_masses.push(molar_mass);
    }
    Ok(molar_masses)
}

/// Build a species-by-element composition matrix. The returned element vector
/// defines the matrix column order and is deterministic (lexicographic).
pub fn create_elem_composition_matrix(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<(DMatrix<f64>, Vec<String>)> {
    let mut hashset_of_elems: BTreeSet<String> = BTreeSet::new();
    let mut vec_of_compositions = Vec::new();
    for formula in vec_of_formulae.iter() {
        // create a unique list of elements from the given formula vector
        let counts = parse_formula(formula.to_string(), groups.clone())?;
        validate_physical_composition(&counts, formula)?;
        vec_of_compositions.push(counts.clone());
        let elements = counts.keys().map(|el| el.clone()).collect::<Vec<_>>();
        hashset_of_elems.extend(elements);
    }
    let unique_vec_of_elems = hashset_of_elems.into_iter().collect::<Vec<_>>();
    let num_rows = unique_vec_of_elems.len();
    let num_cols = vec_of_compositions.len();
    let mut matrix = DMatrix::zeros(num_rows, num_cols);
    for substance_i in 0..vec_of_formulae.len() {
        for j in 0..unique_vec_of_elems.len() {
            let element_j = unique_vec_of_elems[j].clone();
            if let Some(count) = vec_of_compositions[substance_i].get(&element_j) {
                matrix[(j, substance_i)] += *count as f64;
            }
        }
    }
    Ok((matrix.transpose(), unique_vec_of_elems))
}

pub fn create_elem_composition_matrix_and_molar_masses(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<(DMatrix<f64>, Vec<String>, Vec<f64>)> {
    let mut hashset_of_elems: BTreeSet<String> = BTreeSet::new();
    let mut vec_of_compositions = Vec::new();
    let mut vec_of_molar_masses = Vec::new();
    for formula in vec_of_formulae.iter() {
        // create a unique list of elements from the given formula vector
        let (molar_mass, counts) = calculate_molar_mass(formula.to_string(), groups.clone())?;
        vec_of_molar_masses.push(molar_mass);
        vec_of_compositions.push(counts.clone());
        let elements = counts.keys().map(|el| el.clone()).collect::<Vec<_>>();
        hashset_of_elems.extend(elements);
    }
    let unique_vec_of_elems = hashset_of_elems.into_iter().collect::<Vec<_>>();
    let num_rows = unique_vec_of_elems.len();
    let num_cols = vec_of_compositions.len();
    let mut matrix = DMatrix::zeros(num_rows, num_cols);
    for substance_i in 0..vec_of_formulae.len() {
        for j in 0..unique_vec_of_elems.len() {
            let element_j = unique_vec_of_elems[j].clone();
            if let Some(count) = vec_of_compositions[substance_i].get(&element_j) {
                matrix[(j, substance_i)] += *count as f64;
            }
        }
    }
    Ok((matrix.transpose(), unique_vec_of_elems, vec_of_molar_masses))
}

/*
fn main() {
    let formulae = vec!["O12C12", "Na(CLNO)3", "H2O", "H2(CCl3)2", "Na(N2O)2", "C(OOH)2", "C(OOH)2N(ClO)3"];
    for formula in formulae {
        let (molar_mass, element_composition) = calculate_molar_mass(formula.to_string());
        println!("Element counts: {:?}", element_composition);
        println!("Molar mass: {:?} g/mol", molar_mass);
    }


}
*/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_formula() {
        /**/

        let formula = "C6H8O6".to_string();
        let expected_counts = HashMap::from([
            ("C".to_string(), 6),
            ("H".to_string(), 8),
            ("O".to_string(), 6),
        ]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);

        let formula = "Na(NO3)2".to_string();
        let expected_counts = HashMap::from([
            ("Na".to_string(), 1),
            ("N".to_string(), 2),
            ("O".to_string(), 6),
        ]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);

        let formula = "H2O".to_string();
        let expected_counts = HashMap::from([("H".to_string(), 2), ("O".to_string(), 1)]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);

        let formula = "C5H6OOH".to_string();
        let expected_counts = HashMap::from([
            ("C".to_string(), 5),
            ("H".to_string(), 7),
            ("O".to_string(), 2),
        ]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);

        let formula = "Ca(NO3)2".to_string();
        let expected_counts = HashMap::from([
            ("Ca".to_string(), 1),
            ("N".to_string(), 2),
            ("O".to_string(), 6),
        ]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);

        let formula = "C(OOH)2N(ClO)3".to_string();
        let expected_counts = HashMap::from([
            ("C".to_string(), 1),
            ("O".to_string(), 7),
            ("H".to_string(), 2),
            ("N".to_string(), 1),
            ("Cl".to_string(), 3),
        ]);
        assert_eq!(parse_formula(formula, None).unwrap(), expected_counts);
    }

    #[test]
    fn test_calculate_molar_mass() {
        let formula = "H2O(g)".to_string();
        let expected_molar_mass = 18.01528;
        let (molar_mass, _) = calculate_molar_mass(formula, None).unwrap();
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "NaCl".to_string();
        let expected_molar_mass = 58.44;
        let (molar_mass, _) = calculate_molar_mass(formula, None).unwrap();
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "C6H8O6".to_string();
        let expected_molar_mass = 176.12;
        let (molar_mass, _) = calculate_molar_mass(formula, None).unwrap();
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "Ca(NO3)2".to_string();
        let expected_molar_mass = 164.093;
        let (molar_mass, _) = calculate_molar_mass(formula, None).unwrap();
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);
    }

    #[test]
    fn test_calculate_molar_mass_of_vector_of_substances() {
        let vec_of_formulae = vec!["H2O", "NaCl", "C6H8O6", "Ca(NO3)2"];
        let expected_molar_masses = vec![18.01528, 58.44316, 176.12, 164.093];

        let calculated_molar_masses =
            calculate_molar_mass_of_vector_of_subs(vec_of_formulae, None).unwrap();

        for (i, &expected_molar_mass) in expected_molar_masses.iter().enumerate() {
            assert!((calculated_molar_masses[i] - expected_molar_mass).abs() < 1e-2);
        }
    }
    #[test]
    fn test_with_groups() {
        let toluol = "C6H5Me".to_string();
        let expected_counts = HashMap::from([("H".to_string(), 8), ("C".to_string(), 7)]);
        let groups = Some(HashMap::from([(
            "Me".to_string(),
            HashMap::from([("C".to_string(), 1), ("H".to_string(), 3)]),
        )]));
        assert_eq!(parse_formula(toluol, groups).unwrap(), expected_counts);

        let Xylole = "C6H4(Me)2".to_string();
        let expected_counts = HashMap::from([("H".to_string(), 10), ("C".to_string(), 8)]);
        let groups = Some(HashMap::from([(
            "Me".to_string(),
            HashMap::from([("C".to_string(), 1), ("H".to_string(), 3)]),
        )]));
        assert_eq!(parse_formula(Xylole, groups).unwrap(), expected_counts);
    }

    #[test]
    fn test_parse_formula_strips_phase_markers_with_groups() {
        // Phase markers should not interfere with group expansion or stoichiometric parsing.
        let formula = "C6H4(Me)2(g)".to_string();
        let groups = Some(HashMap::from([(
            "Me".to_string(),
            HashMap::from([("C".to_string(), 1), ("H".to_string(), 3)]),
        )]));
        let expected_counts = HashMap::from([("H".to_string(), 10), ("C".to_string(), 8)]);

        assert_eq!(parse_formula(formula, groups).unwrap(), expected_counts);
    }

    #[test]
    fn test_calculate_molar_mass_strips_phase_markers() {
        let formula = "H2O(l)".to_string();
        let expected_molar_mass = 18.01528;
        let (molar_mass, counts) = calculate_molar_mass(formula, None).unwrap();

        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);
        assert_eq!(
            counts,
            HashMap::from([("H".to_string(), 2), ("O".to_string(), 1)])
        );
    }

    #[test]
    fn test_element_matrix() {
        let vec_of_formulae = vec!["H2O", "NaCl", "C3H8", "CH4"]; // 5 elements
        let (matrix, elements) = create_elem_composition_matrix(vec_of_formulae, None).unwrap();
        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 5);
        assert_eq!(elements, vec!["C", "Cl", "H", "Na", "O"]);
        assert_eq!(matrix[(0, 0)], 0.0);
        assert_eq!(matrix[(0, 2)], 2.0);
        assert_eq!(matrix[(0, 4)], 1.0);
        assert_eq!(matrix[(1, 1)], 1.0);
        assert_eq!(matrix[(1, 3)], 1.0);
    }

    #[test]
    fn physical_mass_and_composition_reject_unknown_element_symbols() {
        let counts = parse_formula("Xx2".into(), None).expect("formula tokens remain parseable");
        assert_eq!(counts.get("Xx"), Some(&2));
        assert!(matches!(
            calculate_molar_mass("Xx2".into(), None),
            Err(KineticsError::UnknownElement { symbol, .. }) if symbol == "Xx"
        ));
        assert!(matches!(
            calculate_molar_mass_for_composition(HashMap::from([("Xx".into(), 1)])),
            Err(KineticsError::UnknownElement { symbol, .. }) if symbol == "Xx"
        ));
        assert!(matches!(
            create_elem_composition_matrix(vec!["Xx2"], None),
            Err(KineticsError::UnknownElement { symbol, .. }) if symbol == "Xx"
        ));
    }

    #[test]
    fn empty_formulas_and_zero_atom_counts_are_rejected() {
        for formula in ["", "   ", "(g)", "H0", "Ca(NO3)0"] {
            assert!(
                matches!(
                    parse_formula(formula.to_string(), None),
                    Err(KineticsError::FormulaSyntax { .. })
                ),
                "formula {formula:?} should be rejected"
            );
        }
        assert!(matches!(
            calculate_molar_mass_for_composition(HashMap::new()),
            Err(KineticsError::FormulaSyntax { .. })
        ));
        assert!(matches!(
            calculate_molar_mass_for_composition(HashMap::from([("H".into(), 0)])),
            Err(KineticsError::FormulaSyntax { .. })
        ));
        assert!(matches!(
            parse_formula(
                "Me".into(),
                Some(HashMap::from([(
                    "Me".into(),
                    HashMap::from([("C".into(), 1), ("H".into(), 0)])
                )]))
            ),
            Err(KineticsError::InvalidCustomGroup { .. })
        ));
    }

    #[test]
    fn formula_parser_rejects_non_ascii_before_indexing() {
        assert!(matches!(
            parse_formula("H₂O".into(), None),
            Err(KineticsError::FormulaSyntax { reason, .. }) if reason.contains("ASCII")
        ));
    }

    #[test]
    fn formula_parser_bounds_length_and_nesting_and_checks_count_overflow() {
        let too_long = "H".repeat(MAX_FORMULA_BYTES + 1);
        assert!(matches!(
            parse_formula(too_long, None),
            Err(KineticsError::FormulaLimitExceeded {
                resource: "byte length",
                ..
            })
        ));

        let deep = format!(
            "{}H{}",
            "(".repeat(MAX_FORMULA_NESTING + 1),
            ")".repeat(MAX_FORMULA_NESTING + 1)
        );
        assert!(matches!(
            parse_formula(deep, None),
            Err(KineticsError::FormulaLimitExceeded {
                resource: "nesting depth",
                ..
            })
        ));

        let overflow = format!("H{}H", usize::MAX);
        assert!(matches!(
            parse_formula(overflow, None),
            Err(KineticsError::FormulaSyntax { reason, .. }) if reason.contains("overflow")
        ));
        let multiplier_overflow = format!("(H{})2", usize::MAX);
        assert!(matches!(
            parse_formula(multiplier_overflow, None),
            Err(KineticsError::FormulaSyntax { reason, .. }) if reason.contains("overflow")
        ));
    }

    #[test]
    fn phase_markers_are_removed_only_as_recognized_suffixes() {
        assert_eq!(
            parse_formula("C(OH)2(g)".into(), None).unwrap(),
            HashMap::from([("C".into(), 1), ("O".into(), 2), ("H".into(), 2),])
        );
        assert!(matches!(
            parse_formula("C(g)H".into(), None),
            Err(KineticsError::FormulaSyntax { .. })
        ));
    }

    #[test]
    fn custom_group_expansions_validate_definitions_and_use_longest_prefix_match() {
        let unknown = HashMap::from([("Me".into(), HashMap::from([("Xx".into(), 1)]))]);
        assert!(matches!(
            parse_formula("Me".into(), Some(unknown)),
            Err(KineticsError::UnknownElement { symbol, .. }) if symbol == "Xx"
        ));

        let nested = HashMap::from([
            ("Me".into(), HashMap::from([("Ph".into(), 1)])),
            ("Ph".into(), HashMap::from([("C".into(), 6)])),
        ]);
        assert!(matches!(
            parse_formula("Me".into(), Some(nested)),
            Err(KineticsError::InvalidCustomGroup { reason, .. }) if reason.contains("nested")
        ));

        let prefix_related = HashMap::from([
            ("Me".into(), HashMap::from([("C".into(), 1)])),
            ("Methyl".into(), HashMap::from([("H".into(), 2)])),
        ]);
        assert_eq!(
            parse_formula("Me".into(), Some(prefix_related.clone())).unwrap(),
            HashMap::from([("C".into(), 1)])
        );
        assert_eq!(
            parse_formula("Methyl2".into(), Some(prefix_related)).unwrap(),
            HashMap::from([("H".into(), 4)])
        );
    }

    #[test]
    fn explicit_species_compositions_accept_chemical_identifiers_verbatim() {
        let groups = HashMap::from([
            (
                "C2H5-COOH".into(),
                HashMap::from([("C".into(), 3), ("H".into(), 6), ("O".into(), 2)]),
            ),
            (
                "C2H5-CH2-CH2-CH3".into(),
                HashMap::from([("C".into(), 5), ("H".into(), 12)]),
            ),
        ]);

        assert_eq!(
            parse_formula("C2H5-COOH".into(), Some(groups.clone())).unwrap(),
            HashMap::from([("C".into(), 3), ("H".into(), 6), ("O".into(), 2)])
        );
        assert_eq!(
            parse_formula("C2H5-CH2-CH2-CH3".into(), Some(groups)).unwrap(),
            HashMap::from([("C".into(), 5), ("H".into(), 12)])
        );

        let prefix_named_species = HashMap::from([
            (
                "HMX".into(),
                HashMap::from([
                    ("H".into(), 4),
                    ("N".into(), 8),
                    ("C".into(), 8),
                    ("O".into(), 8),
                ]),
            ),
            (
                "HMXprod".into(),
                HashMap::from([("H".into(), 6), ("C".into(), 1), ("O".into(), 1)]),
            ),
        ]);
        assert_eq!(
            parse_formula("HMX".into(), Some(prefix_named_species.clone())).unwrap(),
            prefix_named_species["HMX"]
        );
        assert_eq!(
            parse_formula("HMXprod".into(), Some(prefix_named_species)).unwrap(),
            HashMap::from([("H".into(), 6), ("C".into(), 1), ("O".into(), 1)])
        );

        let formulas = ["C2H5-COOH", "C2H5-CH2-CH2-CH3"];
        let molar_masses = calculate_molar_mass_of_vector_of_subs(
            formulas.to_vec(),
            Some(HashMap::from([
                (
                    formulas[0].to_string(),
                    HashMap::from([("C".into(), 3), ("H".into(), 6), ("O".into(), 2)]),
                ),
                (
                    formulas[1].to_string(),
                    HashMap::from([("C".into(), 5), ("H".into(), 12)]),
                ),
            ])),
        )
        .unwrap();
        assert!((molar_masses[0] - 74.08).abs() < 0.1);
        assert!((molar_masses[1] - 72.15).abs() < 0.1);
    }
    #[test]
    fn test_hmx() {
        let hmx = HashMap::from([
            ("H".to_string(), 4),
            ("N".to_string(), 8),
            ("C".to_string(), 8),
            ("O".to_string(), 8),
        ]);
        let groups = Some(HashMap::from([("HMX".to_string(), hmx.clone())]));
        let formula = "HMX".to_string();
        let a = parse_formula(formula.clone(), groups.clone()).unwrap();
        println!("{:?}", a);
        assert!(a == hmx);
        let (molar_mass, _) = calculate_molar_mass(formula, groups).unwrap();
        println!("molar mass: {}", molar_mass);
    }
}
