[TOC]

# KiThe

This is a package of structures, functions and databases useful for such areas as chemical thermodynamics, chemical kinetics, as well as modeling of chemical reactors, combustion, processes in shock tubes and rocket engines, propulsion. 

PROJECT NEWS: Kinetics module rewritten in more idiomatic style, more features added 

## Content
- [Kinetics](#usage)
- [Testing](#testing)
- [Contributing](#contributing)
- [To do](#to-do)



## Kinetics
- parse reaction equations into a list of substances 
- parse reaction equations into a stoichiometric matrix, matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the kinetic function,
- calculate of atomic composition, molar masses and matrix of atomic composition.

Let us observe the main structure realizing that features
```rust
pub struct StoichAnalyzer {
    pub reactions: Vec<String>, //a vector of reactions
    pub groups: Option<HashMap<String, HashMap<String, usize>>>, // Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
    pub substances: Vec<String>, // a vector of substances
    pub stecheo_matrx: Vec<Vec<f64>>, //  a vector of vectors of stoichiometric coefficients in each reaction
    pub stecheo_reags: Vec<Vec<f64>>, //  a vector of vectors of stoichiometric coefficients of reactants in each reaction
    pub stecheo_prods: Vec<Vec<f64>>, //  a vector of vectors of stoichiometric coefficients of products in each reaction
    pub G_matrix_reag: Vec<Vec<f64>>,// matrix of powers of concentrations for constructng kinetic equation for forward reaction 
    pub G_matrix_prod: Vec<Vec<f64>>,// matrix of powers of concentrations for constructng kinetic equation for reverse reaction 
    pub matrix_of_elements: Option<DMatrix<f64>>,// matrix of elemental composition
    pub vec_of_molmasses: Option<Vec<f64>>,// vector of molecular masses
}
```
usage
```rust
use KiThe::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
let mut  StoichAnalyzer_instance = StoichAnalyzer::new();
let reactions_: Vec<&str> = vec!["A=2BM)", "B=>A + 3C_DUP", "2B+A=D**0.5"];
let reaction = reactions_.iter().map(|s| s.to_string()).collect();
 StoichAnalyzer_instance.reactions = reaction;
 StoichAnalyzer_instance.search_substances();
 StoichAnalyzer_instance.analyse_reactions();
          
let stecheo_matrx =  StoichAnalyzer_instance.stecheo_matrx;

let result = [
                [-1.0, 2.0, 0.0, 0.0],
                [1.0, -1.0, 3.0, 0.0],
                [-1.0, -2.0, 0.0, 1.0],
            ];
let result: Vec<Vec<f64>> = result.iter().map(|row| row.to_vec()).collect();
assert_eq!(stecheo_matrx, result);
println!("substances: {:?}", StoichAnalyzer_instance.substances);
println!("stecheo_matrx {:?}", stecheo_matrx);
```
 - Calculation of atomic composition, molar masses and matrix of atomic composition. We can use struct StoichAnalyzer
 or just use the more low-level functions
 ```rust
use KiThe::Kinetics::molmass::{calculate_molar_mass, parse_formula, calculate_molar_mass_of_vector_of_subs,
create_elem_composition_matrix};
let formula = "C6H8O6";
let (molar_mass, element_composition) = calculate_molar_mass(formula.to_string(), None); 
 println!("Element counts: {:?}", element_composition);
 println!("Molar mass: {:?} g/mol", molar_mass);
            
  let formula = "Na(NO3)2".to_string();
  let atomic_composition = parse_formula(formula, None);
  println!("{:?}", atomic_composition);


  let vec_of_formulae = vec!["H2O", "NaCl", "C6H8O6", "Ca(NO3)2"];
  let expected_molar_masses = vec![18.01528, 58.44316, 176.12, 164.093];
  let calculated_molar_masses = calculate_molar_mass_of_vector_of_subs(vec_of_formulae, None);
    
for (i, &expected_molar_mass) in expected_molar_masses.iter().enumerate() {
                println!("molar mass: {:?} g/mol", calculated_molar_masses[i]);
                assert!((calculated_molar_masses[i] - expected_molar_mass).abs() < 1e-2);
}

let vec_of_formulae = vec!["H2O", "NaCl", "C3H8", "CH4"]; // 5 elements
let matrix = create_elem_composition_matrix(vec_of_formulae, None);
println!("{}", matrix);
 ```
- crate is equipped with a libraries of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases, so you can
view all libraries and all reactions in every library
of kinetic DB, search reactions by substances and so on. Most important methods below
```rust
use KiThe::Kinetics::kinetics_lib_api::KineticData;
let mut kin_instance = KineticData::new();
            // collecting reaction data for library name lib
let lib = "NUIG";

kin_instance.open_json_files(lib);
            // veiew all reactions in library
kin_instance.print_all_reactions();
      
let reaction1 = kin_instance.search_reactdata_by_reaction_id("1");
println!("reaction1: {:?}", reaction1);
            // search reactions by substances 
kin_instance.search_reaction_by_reagents_and_products((vec!["CO".to_string()])  );
println!("reactions where CO is product: {:?}", kin_instance.FoundReactionsByProducts);
println!("reactions where CO is reagent: {:?}", kin_instance.FoundReactionsByReagents);
```
Kinetic data (Arrhenius parameters, etc) are parsed into stuctures 

```rust
pub struct ReactionData {
    #[serde(rename = "type")]
    reaction_type: ReactionType,
    eq: String,
    pub react: Option<HashMap<String, f64>>,
    #[serde(flatten)]
    data: ReactionKinetics,
}
```
where ReactionType is enum for concrete type of reaction (elementary, fall-off, etc)

-  The module is automatic chemical mechanism constructor and takes as input the name of the library and the vector of substances and then produces the following data:
  1) all reactions of starting substances with each other, and all reactions of all their possible products with each other and with original substances. 
  2) HashMap with kinetic data of all found reactions
```rust
use crate::Kinetics::mechfinder_api::Mechanism_search;
let mut mech_search = Mechanism_search::new(
                vec!["O".to_string(), "NH3".to_string(), "NO".to_string()], // what to search
                "NUIG".to_string(), // library name
            );
    
let (mechanism, reactants, vec_of_reactions) = mech_search.mechfinder_api();
println!("mechanism (reaction ID's) : {:?}", mechanism);
println!("reactants: {:?}", reactants);
println!("reaction data: {:?}", vec_of_reactions);
println!("vector of ReactionData structs with parsed data: {:#?}", mech_search.reactdata);
```
The most general approach is to use stuct KinData whin aggregates the most features of Kinetics module: library parsing,
mechanism construction, structures for different reaction types, stoichiometric calculations, saving and reading data
```rust
pub struct KinData {
    pub shortcut_reactions: Option<Vec<String>>, // vector of reaction shortcut names
    pub map_of_reactions: Option<HashMap<String, Vec<String>>>, // full "address" of reaction {'library':"id of reaction"} group reactions by library names {'library':[reaction_ids in that library]}
    pub vec_of_pairs: Option<Vec<(String, String)>>, // vector of pairs of reaction library names and reaction ids
    pub vec_of_reaction_Values: Option<Vec<Value>>, // vector of reaction data in the form of serde Values
    pub vec_of_reaction_data: Option<Vec<ReactionData>>, // data of all reactions
    pub vec_of_equations: Vec<String>,                        // vector of equations of reactions
    pub substances: Vec<String>,                              // vector of substance names
    pub groups: Option<HashMap<String, HashMap<String, usize>>>, // Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
    pub stecheodata: StoichAnalyzer, // matrix of stoichiometric coefficients and other matrices
}
```
examples of usage 
```rust
use crate::Kinetics::User_reactions::KinData;
    // let our journey begin with a new instance of KinData
let mut kd = KinData::new();
 // set the shortcut reactions for our KineticData instance
// it means we want reactions from Cantera sub-librarie from number 1 to number 10
kd.set_reactions_from_shortcut_range("C1..C10".to_string());
            // searching for reactions in data base
kd.get_reactions_from_shortcuts();
kd.kinetic_main(); // parsing reaction data into structures and stoichometric calculations under one hood
kd.pretty_print_kindata(); // pretty print the reaction data



```
## Testing
Our project is covered by tests and you can run them by standard command
```sh
cargo test
```

## Contributing
If you have any questions, comments or want to contribute, please feel free to contact us at https://github.com/



## To do
- [x] Add libraries of chemical reactions with appropriate methods for processing, searching, and retrieving data.
- [ ] Add libraries of chemical substances...
- [ ] Add numerical methods (may be cpp open source...)




