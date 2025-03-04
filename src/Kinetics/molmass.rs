use nalgebra::DMatrix;
/// Module to calculate the atomic composition and molar mass of a chemical formula
///
///
use std::collections::{HashMap, HashSet};

// Define a struct to hold element data
pub struct Element {
    name: &'static str,
    atomic_mass: f64,
}

// Define a list of elements and their atomic masses
const ELEMENTS: &[Element] = &[
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

    let phases = ["(C)", "(c)", "(L)", "(l)", "(G)", "(g)", "(S)", "(s)", "S)", "s)"];
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
    if let Some(groups) = groups {
        let mut to_remove = Vec::new();

        for (chemical_group, atomic_composition) in groups.iter() {
            // if a group is found in the dictionary
            // we should get rid of it and turn it into regular elements, i.e. Me (methyl) group is converted into {"C":1, "H":3}
            if let Some(&number_of_chemical_groups) = counts.get(chemical_group) {
                to_remove.push(chemical_group.clone());
                for (atom, &quantity) in atomic_composition.iter() {
                    *counts.entry(atom.clone()).or_insert(quantity) +=
                        quantity * number_of_chemical_groups;
                }
            }
        }

        for group in to_remove {
            counts.remove(&group);
        }
    }
    counts
}
fn after_bracket_stoichio(end_bracket: usize, formula: String) -> usize {
    let mut end_of_stoichio_after_bracket = end_bracket + 1;
    while end_of_stoichio_after_bracket < formula.len()
        && formula
            .chars()
            .nth(end_of_stoichio_after_bracket)
            .unwrap()
            .is_digit(10)
    {
        end_of_stoichio_after_bracket += 1;
    }
    let stoichio: usize = formula[end_bracket..end_of_stoichio_after_bracket]
        .parse()
        .unwrap();
    stoichio
}

// Function to parse a chemical formula and return a HashMap of elements and their counts. Argument groups are optional. It is
// needed if formula contains shecial names for chemical groups like Me, Ph, etc In that case this argument should contain the names of these groups
// and there atomic composition { "Me":{"C":1, "H":3}}
pub fn parse_formula(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> HashMap<String, usize> {
    let mut counts = HashMap::new();
    let mut i = 0;
    let mut formula = formula.replace(" ", "");
    let initial_formula = formula.clone();
    println!("PARSING FORMULA: {}", formula);
    formula = filter_phases_marks(&formula);
    let mut start_bracket = 0;
    let mut end_bracket = 0;
    while i < formula.len() {
        let start = i;
        println!("new cicle with position {}", start);
        // proceedeing cases with brackets in formula
        // let's find positions of element inside brackets
        let mut j = i;
        if formula.chars().nth(i).unwrap() == '(' {
            println!("opening brecket found at position {}, deleing brackets", i);
            // move bracket out

            formula.replace_range(i..i + 1, "");
            // position of opening bracket
            start_bracket = j;
            j += 1;
            // move to position of closing bracket
            while i < formula.len() && j < formula.len() && formula.chars().nth(j).unwrap() != ')' {
                j += 1;
            }
            // position of closing bracket
            end_bracket = j;
            println!("closing brecket found at position {}, deleing brackets", j);
            // move bracket out

            formula.replace_range(j..j + 1, "");
            println!(
                "length of formula after deleting brackets {}",
                formula.len()
            );
        }
        //  println!("start bracket position {}, end bracket position {}", start_bracket, end_bracket);
        // vector of element names
        let mut element_names = ELEMENTS
            .iter()
            .map(|element| element.name)
            .collect::<Vec<_>>();
        if let Some(groups) = groups.as_ref() {
            let add_groups_to_elements: Vec<&str> =
                groups.keys().map(|group| group.as_str()).collect();
            element_names.extend(add_groups_to_elements);
        }

        // if letter in the current position is uppercase it is an element or a first latter of an element name
        if formula.chars().nth(i).unwrap().is_uppercase() && i < formula.len() {
            println!(
                "Find uppercase letter: {},  move to the next position {}",
                formula.chars().nth(i).unwrap(),
                i + 1
            );
            i += 1;

            // let's find out is it a two-latter element name and if it is, move to the next position
            // error handling in case this is last element symbol in formula like O in formula H2O
            if let Some(c) = formula.chars().nth(i) {
                if c.is_lowercase() {
                    println!(
                        "Find lowercase letter: {}, move to the next position {}",
                        formula.chars().nth(i).unwrap(),
                        i + 1
                    );
                    i += 1;
                } else if c.is_uppercase() && !element_names.contains(&c.to_string().as_str()) {
                    // If current position is a upperrcase and it is a second latter of an element name
                    // we must turn it into lowercase and move to the next position
                    formula.replace_range(
                        i..i + 1,
                        formula
                            .chars()
                            .nth(i)
                            .unwrap()
                            .to_lowercase()
                            .to_string()
                            .as_str(),
                    );
                    i += 1;
                }
            } else {
                println!("Index out of bounds: {}", i);
            } // end of error handling
        }

        // get an element name
        #[allow(unused_variables)]
        let element: &str;

        if start < i {
            element = &formula[start..i];
        } else {
            println!("break");
            break;
        }

        //let element = &formula[start..i];
        println!("element found: {}", element);
        // if element is empty we should panic
        if element.is_empty() {
            panic!("element is empty");
        }

        let mut count = 1;

        // if we have no elements in brackets or current element is out of brackets
        if (start_bracket == 0 && end_bracket == 0 && i < formula.len())
            || i < start_bracket
            || i > end_bracket
        {
            println!("element not in brackets {}", element);
            if let Some(n) = formula.chars().nth(i) {
                if n.is_digit(10) {
                    let mut end = i + 1;
                    while end < formula.len() && formula.chars().nth(end).unwrap().is_digit(10) {
                        println!("end: {}", end);
                        end += 1;
                    }
                    count = formula[i..end].parse().unwrap();
                    i = end;
                } // end of if let n
                println!("stoichiometric {} for element: {}", count, element);
            } // end of if  start_bracket==0 && end_bracket==0 
        }
        // end of if start_bracket==0 && end_bracket==0
        // if we have elements in brackets and current element is in brackets
        else if start_bracket != 0
            && end_bracket != 0
            && i < formula.len()
            && i > start_bracket
            && i <= end_bracket
        {
            let i_is_digit = formula.chars().nth(i).unwrap().is_digit(10);
            // if element is not a stoichomical coefficient it is an element
            if !i_is_digit {
                println!("it is a digit?: {}", i_is_digit);
                println!(
                    "bracket positions {}, {}, current position: {}",
                    start_bracket, end_bracket, i
                );
                println!("element {} inside brackets", element);
                let stoichio: usize = after_bracket_stoichio(end_bracket, formula.clone());
                count = stoichio;

                println!("new stoichio: {}, for element: {}", count, element);
            // found stoichio inside brackets
            } else if i_is_digit {
                // in sunstances like (A3..B2)2 stoichiometric coefficient = coefficient inside brackets * stoichiometric coefficient outside brackets
                println!(
                    "find digit: {} inside brackets",
                    formula.chars().nth(i).unwrap()
                );
                println!("it is a digit?: {}", i_is_digit);
                println!(
                    "bracket positions {}, {}, current position: {}",
                    start_bracket, end_bracket, i
                );
                println!("element {} inside brackets", element);
                let stoichio_after_bracket: usize =
                    after_bracket_stoichio(end_bracket, formula.clone());
                println!("stoichiometric after bracket: {}", stoichio_after_bracket);
                let mut end = i + 1;
                // condition end<end_bracket added for cases like (..H6)2 to avoid parsing then into H62
                while end < formula.len()
                    && formula.chars().nth(end).unwrap().is_digit(10)
                    && end < end_bracket
                {
                    println!("end: {}", end);
                    end += 1;
                }
                // condition if i < end_bracket added for cases like (..H)2 to avoid counting 2 twice as stoichio_after_elemen and as stoichio_after_bracket
                let stoichio_after_element: usize = if i < end_bracket {
                    formula[i..end].parse().unwrap()
                } else {
                    1
                };
                println!("stoichiometric after element: {}", stoichio_after_element);
                i = end;

                count = stoichio_after_bracket * stoichio_after_element;

                println!("count: {}, i {}", count, i);
                //step forward
                //
            } // end of if i_is_digit
        } // end of if start_bracket!=0 && end_bracket!=0 && i < formula.len() && i > start_bracket && i <= end_bracket

        // update hashmap
        *counts.entry(element.to_string()).or_insert(0) += count;
        println!(
            "hashmap of elements updated: {:?} added {} atoms {} at position {}",
            &counts, count, element, i
        );
        if i >= formula.len() {
            println!("end of parsing formula {}", &initial_formula);
            break;
        } else {
            println!(
                "position: {} while length of formula: {} - go to next loop",
                i,
                formula.len()
            );
        }
    } //  end of while loop
    let counts = handle_groups(counts, groups);
    counts
} // end of parse_formula

// Function to calculate the molar mass of a substance given its chemical formula
pub fn calculate_molar_mass(
    formula: String,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> (f64, HashMap<String, usize>) {
    let counts = parse_formula(formula, groups);
    let mut molar_mass = 0.0;
    for (element, count) in counts.clone() {
        for e in ELEMENTS {
            if e.name == element {
                println!("found element: {}, number of atoms  {}", e.name, count);
                molar_mass += e.atomic_mass * count as f64;
                break;
            }
        }
    }
    (molar_mass, counts)
}

// Function to calculate the molar mass of a vector of chemical formulas
pub fn calculate_molar_mass_of_vector_of_subs(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) -> Vec<f64> {
    println!("\n___________CALCULATE MOLAR MASS OF VECTOR OF SUBS___________");
    let mut molar_masses = Vec::new();
    for formula in vec_of_formulae {
        let counts = parse_formula(formula.to_string(), groups.clone());
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
    println!("___________CALCULATE MOLAR MASS OF VECTOR OF SUBS ENDED___________");
    molar_masses
}

pub fn create_elem_composition_matrix(
    vec_of_formulae: Vec<&str>,
    groups: Option<HashMap<String, HashMap<String, usize>>>,
) ->( DMatrix<f64>, Vec<String>) {
    println!("\n___________CREATE ELEMENTS COMPOSITION MATRIX___________");
    let mut hashset_of_elems: HashSet<String> = HashSet::new();
    let mut vec_of_compositions = Vec::new();
    for formula in vec_of_formulae.iter() {
        // create a unique list of elements from the given formula vector
        let counts = parse_formula(formula.to_string(), groups.clone());
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
    println!("\n___________CREATE ELEMENTS COMPOSITION MATRIX ENDED___________");
    (matrix.transpose(), unique_vec_of_elems)
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
        assert_eq!(parse_formula(formula, None), expected_counts);

        let formula = "Na(NO3)2".to_string();
        let expected_counts = HashMap::from([
            ("Na".to_string(), 1),
            ("N".to_string(), 2),
            ("O".to_string(), 6),
        ]);
        assert_eq!(parse_formula(formula, None), expected_counts);

        let formula = "H2O".to_string();
        let expected_counts = HashMap::from([("H".to_string(), 2), ("O".to_string(), 1)]);
        assert_eq!(parse_formula(formula, None), expected_counts);

        let formula = "C5H6OOH".to_string();
        let expected_counts = HashMap::from([
            ("C".to_string(), 5),
            ("H".to_string(), 7),
            ("O".to_string(), 2),
        ]);
        assert_eq!(parse_formula(formula, None), expected_counts);
    }

    #[test]
    fn test_calculate_molar_mass() {
        let formula = "H2O(g)".to_string();
        let expected_molar_mass = 18.01528;
        let (molar_mass, _) = calculate_molar_mass(formula, None);
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "NaCl".to_string();
        let expected_molar_mass = 58.44;
        let (molar_mass, _) = calculate_molar_mass(formula, None);
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "C6H8O6".to_string();
        let expected_molar_mass = 176.12;
        let (molar_mass, _) = calculate_molar_mass(formula, None);
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);

        let formula = "Ca(NO3)2".to_string();
        let expected_molar_mass = 164.093;
        let (molar_mass, _) = calculate_molar_mass(formula, None);
        assert!((molar_mass - expected_molar_mass).abs() < 1e-2);
    }

    #[test]
    fn test_calculate_molar_mass_of_vector_of_substances() {
        let vec_of_formulae = vec!["H2O", "NaCl", "C6H8O6", "Ca(NO3)2"];
        let expected_molar_masses = vec![18.01528, 58.44316, 176.12, 164.093];

        let calculated_molar_masses = calculate_molar_mass_of_vector_of_subs(vec_of_formulae, None);

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
        assert_eq!(parse_formula(toluol, groups), expected_counts);

        let Xylole = "C6H4(Me)2".to_string();
        let expected_counts = HashMap::from([("H".to_string(), 10), ("C".to_string(), 8)]);
        let groups = Some(HashMap::from([(
            "Me".to_string(),
            HashMap::from([("C".to_string(), 1), ("H".to_string(), 3)]),
        )]));
        assert_eq!(parse_formula(Xylole, groups), expected_counts);
    }
    #[test]
    fn test_element_matrix() {
        let vec_of_formulae = vec!["H2O", "NaCl", "C3H8", "CH4"]; // 5 elements
        let matrix = create_elem_composition_matrix(vec_of_formulae, None).0;
        println!("{}", matrix);
        assert_eq!(matrix.nrows(), 4);
        assert_eq!(matrix.ncols(), 5);
    }
}
