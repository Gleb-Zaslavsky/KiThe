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
use std::collections::{HashMap, HashSet};

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

fn filter_phases_marks(formula: &str) -> String {
    let mut formula = formula.to_string();

    let phases = [
        "(C)", "(c)", "(L)", "(l)", "(G)", "(g)", "(S)", "(s)", "S)", "s)",
    ];
    for phase in phases {
        formula = formula.replace(phase, "");
    }
    formula
}
// Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
// so we need to convert them into regular elements
fn handle_groups(
    mut counts: HashMap<String, usize>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> HashMap<String, usize> {
    //println!("DEBUG handle_groups: Input counts: {:?}", counts);
    //println!("DEBUG handle_groups: Groups: {:?}", groups);

    if let Some(groups) = groups {
        let mut to_remove = Vec::new();

        for (chemical_group, atomic_composition) in groups.iter() {
            // if a group is found in the dictionary
            // we should get rid of it and turn it into regular elements, i.e. Me (methyl) group is converted into {"C":1, "H":3}
            if let Some(&number_of_chemical_groups) = counts.get(chemical_group) {
                // println!("DEBUG handle_groups: Found group '{}' with count {}, composition: {:?}",
                //  chemical_group, number_of_chemical_groups, atomic_composition);
                to_remove.push(chemical_group.clone());
                for (atom, &quantity) in atomic_composition.iter() {
                    // let old_count = counts.get(atom).unwrap_or(&0);
                    // let new_count = old_count + quantity * number_of_chemical_groups;
                    // println!("DEBUG handle_groups: Adding {} {} atoms (was {}, now {})",
                    //  quantity * number_of_chemical_groups, atom, old_count, new_count);
                    *counts.entry(atom.clone()).or_insert(0) +=
                        quantity * number_of_chemical_groups;
                }
            }
        }

        for group in to_remove {
            //   println!("DEBUG handle_groups: Removing group '{}'", group);
            counts.remove(&group);
        }
    }

    //  println!("DEBUG handle_groups: Final counts: {:?}", counts);
    counts
}
fn formula_error(formula: &str, message: impl Into<String>) -> KineticsError {
    KineticsError::InvalidReactionData(format!("invalid formula `{}`: {}", formula, message.into()))
}

fn char_at(formula: &str, idx: usize) -> KineticsResult<char> {
    formula
        .chars()
        .nth(idx)
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

fn after_bracket_stoichio(end_bracket: usize, formula: String) -> KineticsResult<usize> {
    let mut end_of_stoichio_after_bracket = end_bracket + 1;
    while end_of_stoichio_after_bracket < formula.len()
        && char_at(&formula, end_of_stoichio_after_bracket)?.is_digit(10)
    {
        end_of_stoichio_after_bracket += 1;
    }
    parse_usize_slice(&formula, end_bracket, end_of_stoichio_after_bracket)
}

/// Parse a numeric suffix that starts at `start` and return both the value and the first
/// index after the suffix.
fn parse_stoichio_suffix(formula: &str, start: usize) -> KineticsResult<(usize, usize)> {
    let mut end = start + 1;
    while end < formula.len() && char_at(formula, end)?.is_digit(10) {
        end += 1;
    }
    Ok((parse_usize_slice(formula, start, end)?, end))
}

fn merge_counts(
    target: &mut HashMap<String, usize>,
    source: HashMap<String, usize>,
    multiplier: usize,
) {
    for (atom, count) in source {
        *target.entry(atom).or_insert(0) += count * multiplier;
    }
}

/// Parse one formula segment until a closing parenthesis or the end of the string.
fn parse_formula_segment(
    formula: &str,
    mut i: usize,
    group_names: &[String],
    in_brackets: bool,
    initial_formula: &str,
) -> KineticsResult<(HashMap<String, usize>, usize)> {
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
            let (inner_counts, next_i) =
                parse_formula_segment(formula, i + 1, group_names, true, initial_formula)?;
            i = next_i;

            let multiplier = if i < formula.len() && char_at(formula, i)?.is_digit(10) {
                let (multiplier, next_after_multiplier) = parse_stoichio_suffix(formula, i)?;
                i = next_after_multiplier;
                multiplier
            } else {
                1
            };

            merge_counts(&mut counts, inner_counts, multiplier);
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

            *counts.entry(group_name.clone()).or_insert(0) += multiplier;
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

            *counts.entry(element.to_string()).or_insert(0) += multiplier;
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

// Function to parse a chemical formula and return a HashMap of elements and their counts. Argument groups are optional. It is
// needed if formula contains shecial names for chemical groups like Me, Ph, etc In that case this argument should contain the names of these groups
// and there atomic composition { "Me":{"C":1, "H":3}}
pub fn parse_formula(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<HashMap<String, usize>> {
    let mut formula = formula.replace(" ", "");
    let initial_formula = formula.clone();
    formula = filter_phases_marks(&formula);

    let mut group_names = groups
        .as_ref()
        .map(|group_map| {
            let mut names = group_map.keys().cloned().collect::<Vec<_>>();
            names.sort_by(|a, b| b.len().cmp(&a.len()));
            names
        })
        .unwrap_or_default();

    // Element symbols remain available for the parser through plain uppercase/lowercase matching,
    // while custom groups are matched first so longer symbolic names do not get split apart.
    if group_names.is_empty() {
        group_names = Vec::new();
    }

    let (counts, _) = parse_formula_segment(&formula, 0, &group_names, false, &initial_formula)?;
    Ok(handle_groups(counts, groups))
} // end of parse_formula

// Function to calculate the molar mass of a substance given its chemical formula
pub fn calculate_molar_mass(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<(f64, HashMap<String, usize>)> {
    let counts = parse_formula(formula, groups)?;
    let mut molar_mass = 0.0;
    for (element, count) in counts.clone() {
        for e in ELEMENTS {
            if e.name == element {
                molar_mass += e.atomic_mass * count as f64;
                break;
            }
        }
    }
    Ok((molar_mass, counts))
}

pub fn calculate_molar_mass_for_composition(counts: HashMap<String, usize>) -> f64 {
    let mut molar_mass = 0.0;
    for (element, count) in counts.clone() {
        for e in ELEMENTS {
            if e.name == element {
                molar_mass += e.atomic_mass * count as f64;
                break;
            }
        }
    }
    molar_mass
}

// Function to calculate the molar mass of a vector of chemical formulas
pub fn calculate_molar_mass_of_vector_of_subs(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<Vec<f64>> {
    let mut molar_masses = Vec::new();
    for formula in vec_of_formulae.iter() {
        let counts = parse_formula(formula.to_string(), groups.clone())?;
        let mut molar_mass = 0.0;
        for (element, count) in counts {
            for e in ELEMENTS {
                if e.name == element {
                    molar_mass += e.atomic_mass * count as f64;
                    break;
                }
            }
        }
        molar_masses.push(molar_mass);
    }
    Ok(molar_masses)
}

pub fn create_elem_composition_matrix(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> KineticsResult<(DMatrix<f64>, Vec<String>)> {
    let mut hashset_of_elems: HashSet<String> = HashSet::new();
    let mut vec_of_compositions = Vec::new();
    for formula in vec_of_formulae.iter() {
        // create a unique list of elements from the given formula vector
        let counts = parse_formula(formula.to_string(), groups.clone())?;
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
    let mut hashset_of_elems: HashSet<String> = HashSet::new();
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
        let matrix = create_elem_composition_matrix(vec_of_formulae, None)
            .unwrap()
            .0;
        println!("{}", matrix);
        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 5);
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
