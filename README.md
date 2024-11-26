[TOC]

# KiThe

This is a package of structures, functions and databases useful for such areas as chemical thermodynamics, chemical kinetics, as well as modeling of chemical reactors, combustion, processes in shock tubes and rocket engines, propulsion. 

## Content
- [Usage](#usage)
- [Testing](#testing)
- [Contributing](#contributing)
- [To do](#to-do)



## Usage
- parse reaction equations into a list of substances 
```rust
use KiThe::reaction_analyzer::ReactionAnalyzer;
let mut ReactionAnalyzer_instance = ReactionAnalyzer::new();
let reactions_: Vec<&str> = vec!["A=B", "B->A + 3C", "2B+A=D"];
let reaction = reactions_.iter().map(|s| s.to_string()).collect();
ReactionAnalyzer_instance.reactions = reaction;
ReactionAnalyzer_instance.search_substances();
println!("substances: {:?}", ReactionAnalyzer_instance.substances);
```
- parse reaction equations into a stoichiometric matrix, matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the 
 kinetic function,
```rust
use KiThe::reaction_analyzer::ReactionAnalyzer;
let mut ReactionAnalyzer_instance = ReactionAnalyzer::new();
let reactions_: Vec<&str> = vec!["A=B", "B->A + 3C", "2B+A=D"];
let reaction = reactions_.iter().map(|s| s.to_string()).collect();
ReactionAnalyzer_instance.reactions = reaction;
ReactionAnalyzer_instance.search_substances();
ReactionAnalyzer_instance.analyse_reactions();
println!("substances: {:?}", ReactionAnalyzer_instance.substances);
println!("{:?}", ReactionAnalyzer_instance);
```
- crate is equipped with a libraries of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases, so you can
view all libraries and all reactions in every library
of kinetic DB, search reactions by substances and so on. Most important methods below
```rust
let mut kin_instance = KineticData::new();
// collecting reaction data for library name lib
kin_instance.open_json_files(lib);
// veiew all reactions in library
kin_instance.print_all_reactions();
// returns reaction ID and reaction data (parsed from json) for given reaction equation
kin_instance.search_reaction_by_equation(equation)
// search reactions by substances 
kin_instance.search_reaction_by_reagents_and_products(reagents)
```
-  The module is automatic chemical mechanism constructor and takes as input the name of the library and the vector of substances and then produces the following data:
  1) all reactions of starting substances with each other, and all reactions of all their possible products with each other and with original substances. 
  2) HashMap with kinetic data of all found reactions
```rust
    let mut mech_search = Mechanism_search::new(
            vec!["O".to_string(), "NH3".to_string(), "NO".to_string()],
            "NUIG".to_string(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );

        let (mechanism, reactants, vec_of_reactions) = mech_search.mechfinder_api();
```
 - Calculation of atomic composition and molar mass
 ```rust
use KiThe::molmass::calculate_molar_mass;
 let formula = "C6H8O6";
let (molar_mass, element_composition) = calculate_molar_mass(formula.to_string()); 
 println!("Element counts: {:?}", element_composition);
 println!("Molar mass: {:?} g/mol", molar_mass);
use KiThe::molmass::parse_formula;
let formula = "Na(NO3)2".to_string();
let atomic_composition = parse_formula(formula);
println!("{:?}", atomic_composition);
 ```
## Testing
Our project is covered by tests and you can run them by standard command
```sh
cargo test
```

## Contributing
If you have any questions, comments or want to contribute, please feel free to contact us at https://github.com/



## To do
here's a roadmap for further development of the crate 
1) Add libraries of thermochemical properties (S, H, Cp), as well as transport properties (viscosity, thermal conductivity, diffusion) and modules for their processing (we have all this as Python code and should be rewritten in Rust). The modules should support different formats of heat capacity polynomials, namely NASA and NIST formats. It is desirable to add a parser of the online database NIST, which in case of absence of the required substance in the local database would look for data on it in NIST (all this we have in the form of code in Python and should be rewritten in Rust).
2) Chemical equilibrium problem. The initial composition of substances and a pair of thermodynamic parameters (usually P, V or H, P) are given and the equilibrium composition and other thermodynamic functions must be found. Such an algorithm is implemented in some commercial packages such as Terra, HSC, etc. and is well described in the works of Trusov and his students. 
3) Add different types of Arrhenius functions (‘fall-off’, ‘three-body’, etc) presented in kinetic libraries. And kinetic functions should be both in the form of regular Rust functions and in the form of symbolic RustedSciThe functions for further calculation of symbolic jacobians.
4) Add a constructor of kinetic expressions: we take the Arrhenius function from point (3) of this roadmap, and the reaction_analyzer module will give us the degrees of kinetic polynomials and stoichiometric coefficients
5) At this stage, we have everything we need to solve the simplest computational problem of a ‘0-dimensional’ reactor: an ideal mixing reactor maintained at constant temperature - "pure kinetics" no mass heat transfer. It is necessary to provide the output of the calculation result in the form of tables, detailed log of the calculation progress and plots.
6) Then it is logical to proceed to the ‘0-dimensional’ problem with initial condition for the ideal mixing reactor with variable temperature, here the heat balance equation is solved, so this stage will rely on the results of all previous stages. The equation systems of steps 4 and 5 will utilise one of the solvers for the IVP from the RST crate. The model solved in this step is also used in the modelling of shock tubes.
7)  Stationary one-dimensional problem with initial condition corresponds to the steady-state regime in a tubular reactor, combustion of condensed phase in models of solid rocket engines, etc. The type of equations is more complicated because diffusion and heat conduction are involved. It is solved by one of the algorithms for solving stiff and moderately stiff equations (BDF, Radau, Backward Euler).
8) Stationary one-dimensional problem with boundary conditions corresponds to the steady-state regime in a tubular reactor, combustion of the gas phase in models of solid rocket engines, etc. The type of equations is more complicated because diffusion and heat conduction are involved. It is solved by one of the types of modified Newton's method or by one of the variants of the collocation method.

solving these problems will cover up to 80% of the computational needs of many kineticists and chemical modellers...
