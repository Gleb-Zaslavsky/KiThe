[TOC]

# KiThe

This is a package of structures, functions and databases useful for such areas as chemical thermodynamics, chemical kinetics, as well as modeling of chemical reactors, combustion, processes in shock tubes and rocket engines, propulsion. 

PROJECT NEWS: BVP solver for conbustion/plug-flow steady state gas-phase reactor now called via CLI menu
## Content
- [Kinetics](#Kinetics)
- [Thermodynamics](#Thermodynamics)
- [Chemical thermodynamics](#Chemical_thermodynamics)
- [NIST scrapper](#NIST_scrapper)
- [Testing](#Testing)
- [Contributing](#contributing)
- [To do](#to-do)

## Features
When install as library:
* Chemical kinetics
    * crate is equipped with a libraries of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases;
    * searching inside local kinetic libraries by reagents, products, etc. 
    * parsing reaction equations into a list of substances; 
    * parsing reaction equations into a stoichiometric matrix, matrix of coefficients of direct reactions and matrix  of coefficients of reverse reactions, matrix of degrees of concentration for the kinetic function;
    * calculation of atomic composition, molar masses and matrix of atomic composition;
    * Automatic (for found in libs reactions) and by function call constucting of kinetic functions of all main types (elemntary, fall-off, troe, etc) - both rust functions and symbolic expressions;
    * Automatic chemical mechanism constructor (say you have some reagents and want to find all possible reactions 
    between original species and all their products)
* Thermodynamics
    * many libraries on board with thermodynamics and transport properties (NASA, NASA-CEA, NIST, Aramco transport, etc.) and handlers for them;
    * search substances thermodynamics and heat-mass transfer data through all libraries with storing of data in a structure;
    * Calculaton of Gibbs free energy of a given mixure of substances (numerical result at given T, P, concentration) and symbolic Gibbs free energy of a given mixure of substances (symbolic result at given T, P, concentration).
    * automatic NIST parser for thermochemical data;
* Gaseos combustion/Plug-flow steady-state 1D boundary value problem.
* IVP problem for multiple solid state kinetic models with constant or linear increasing temperature

 When install as executable:
  CLI instrument to solve some of the heat-mass transfer, combustion, chemical engeneering problems. Now available:
  * BVP solver for gas phase steady-state combustion/plug-flow with constant mass velocity BVP problem   
  * IVP problem for multiple solid state kinetic models with constant or linear increasing temperature

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
  let (matrix, vec_of_formulae) = create_elem_composition_matrix(vec_of_formulae, None);
    
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
println!("reaction eq's: {:?}", vec_of_reactions);
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
mechsnism constracting
```rust
use KiThe::Kinetics::User_reactions::KinData;

let mut kd = KinData::new();
let task_substances = vec!["O".to_string(), "NH3".to_string(), "NO".to_string()];
let task_library = "NUIG".to_string();
kd.construct_mechanism(task_substances, task_library);
kd.kinetic_main(); // parsing reaction data into structures and stoichometric calculations under one hood
kd.pretty_print_kindata();
println!("vector of reactions \n\n {:#?}", kd.vec_of_equations);
println!("vector of substances \n\n {:#?}", kd.substances);

```
## Thermodynamics
search data: Cp, dH, dS, thermal conductivity, viscosity, diffusion coefficient, etc
```rust
use KiThe::Thermodynamics::thermo_lib_api::ThermoData;
let mut thermo_data = ThermoData::new();
println!("Libraries on board {:?} \n \n", thermo_data.AllLibraries);
let this_lib_subs = thermo_data.subs_of_this_lib("Cantera_nasa_base_gas");
 println!("subs of this lib {:?} \n \n", this_lib_subs);
 // you may want to search only in chosen libs
let allowed_libs: Vec<String> = vec!["Aramco_transport".to_string(), "CEA".to_string()];
                          // let us find out what info we have about CO
thermo_data.search_libs_for_subs(vec!["CO".to_string()], None);
println!(
    "hashmap_of_thermo_data {:?}",
    thermo_data.hashmap_of_thermo_data
        );
thermo_data.pretty_print_thermo_data();

```
NASA-CEA lib contains info for calculating thermal dependence of thermal conductivity and viscosity. 
CEAdata module contains corresponding structures and methods for calculating this values.
```rust
use KiThe::Thermodynamics::DBhandlers::CEAdata::CEAdata;
let thermo_data = ThermoData::new();
let sublib = thermo_data.LibThermoData.get("CEA").unwrap();
let CO_data = sublib.get("CO").unwrap();
let mut CEA = CEAdata::new();
CEA.from_serde(CO_data.clone()).unwrap();
CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();
CEA.parse_coefficients().unwrap();
CEA.extract_coefficients(500.0).unwrap();
let lambda = CEA.calculate_Lambda(500.0).unwrap();
let visc = CEA.calculate_Visc(500.0).unwrap();
println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
```
Aramco_transport lib contains info for calculating thermal dependence of thermal conductivity, viscosity and diffusion coefficient. 
TRANSPORTdata module contains corresponding structures and methods for calculating this values.
```rust
use  KiThe::Thermodynamics::DBhandlers::TRANSPORTdata::TransportData;
use KiThe::Thermodynamics::DBhandlers::NASAdata::NASAdata;
let thermo_data = ThermoData::new();
let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
let CO_data = sublib.get("CO").unwrap();
println!("CO_data: {}", CO_data);
// creating instance of the transport TransportData structure 
let mut tr = TransportData ::new();
// filling transport data from serde json
tr.from_serde(CO_data.clone());
// setting units
tr.set_M_unit(Some("g/mol".to_owned())  );
tr.set_P_unit(Some("atm".to_owned()));
tr.set_V_unit( Some("mkPa*s".to_owned()));
tr. set_lambda_unit(Some("mW/m/K".to_owned()));
let T = 473.15; // K 
tr.P = 1.0;
tr.M = 28.0; // g/mol
// calculate viscosity
tr.calculate_Visc(T);
assert_relative_eq!( tr.V, 25.2, epsilon = 5.0);
println!("Viscosity: {:?} mkPa*s", tr.V);
// we need heat capacity - so we shall take it from another library     
let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
let CO_data = sublib.get("CO").unwrap();
let mut NASA = NASAdata::new();
NASA.from_serde(CO_data.clone());
NASA.extract_coefficients(T);
NASA.calculate_Cp_dH_dS(T);
let Cp = NASA.Cp;
println!("Cp: {}", Cp, );
// for density calculation we use now ideal gas equations
 let ro = (tr.P*101325.0)*(tr.M/1000.0)/(R*T);
 // thermal conductivity
let L = tr.calculate_Lambda(Cp, ro, T);
println!("Lambda: {}",L);
```
NASA _gas lib contains info for calculating thermal dependence of heat capacity, enthalpy and entropy.
```rust
use KiThe::Thermodynamics::DBhandlers::NASAdata::NASAdata;
let thermo_data = ThermoData::new();
let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
let CO_data = sublib.get("CO").unwrap();
let mut NASA = NASAdata::new();

NASA.calculate_Cp_dH_dS(400.0);
let Cp = NASA.Cp;
let dh = NASA.dh;
let ds = NASA.ds;

println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
        
let t = 400.0;
NASA.create_closures_Cp_dH_dS();

let Cp_fun = &NASA.C_fun;
let dh_fun = &NASA.dh_fun;
let ds_fun = &NASA.ds_fun;
assert_relative_eq!((Cp_fun)(t), NASA.Cp, epsilon = 1e-6);
assert_relative_eq!((dh_fun)(t), NASA.dh, epsilon = 1e-6);
assert_relative_eq!((ds_fun)(t), NASA.ds, epsilon = 1e-6);
// symbolic expressions
NASA.create_sym_Cp_dH_dS();
let Cp_sym = &NASA.Cp_sym;
let Cp_T = Cp_sym.lambdify1D();
let Cp_value = Cp_T(400.0);
assert_relative_eq!(Cp_value, NASA.Cp, epsilon = 1e-6);
let dh_sym = &NASA.dh_sym;
let dh_T = dh_sym.lambdify1D();
let dh_value = dh_T(400.0);
assert_relative_eq!(dh_value, NASA.dh, epsilon = 1e-6);
        let ds_sym = &NASA.ds_sym;
        let ds_T = ds_sym.lambdify1D();
        let ds_value = ds_T(400.0);
        assert_relative_eq!(ds_value, NASA.ds, epsilon = 1e-6);
```
## Chemical_thermodynamics
- calculating Gibbs free energy for gasess and solids (and their mixtures). Only ideal gases are supported now.

```rust
 // calculating Gibbs free energy withoun concentration correction RT*ln(P/P°) + RT*ln(w_i)
let subs =  vec!["CO".to_string(), "CO2".to_string()];
// calling instance of strucutre SubsData created to search substances data in the databases, store 
// search results and calculate thermo properties
let mut subdata = SubsData::new();
// Set up library priorities
subdata.set_multiple_library_priorities(
        vec!["NASA_gas".to_string()],
        LibraryPriority::Priority,);
subdata.substances = subs.clone();
    // Perform the search
subdata.search_substances();
    // Calling instance of structure Thermodynamics to calculate thermo dG 
let mut thermo = Thermodynamics::new();  
    // Set up basic parameters
thermo.set_T(400.0);
thermo.set_P(101325.0, None); 
let concentration = Some(vec![0.5, 0.5]);
    // Add test substances
thermo.vec_of_substances = subs.clone();   
    // Set up phases
thermo.vec_of_phases.insert(subs[0].clone(), Some(Phases::Gas));
thermo.vec_of_phases.insert(subs[1].clone(), Some(Phases::Gas));
    // savng  the search results in the structure Thermodynamics
thermo.subdata = subdata;
    // Calculate Gibbs free energy
thermo.calculate_Gibbs_free_energy(400.0,  concentration.clone());
    // Calculate symbolic Gibbs free energy
thermo.calculate_Gibbs_sym(400.0, Some(vec![Expr::Var("w1".to_string()), Expr::Var("w2".to_string())]));
thermo.set_P_to_sym();
// getting the results of the calculation
thermo.calculate_Gibbs_fun( 400.0);
let map_of_gibbs = thermo.dG.clone();
let map_of_gibbs_sym = thermo.dG_sym.clone();
let map_of_gibbs_fun = thermo.clone().dG_fun;
for substance in &thermo.vec_of_substances {
            println!("substance: {:?}", substance);
            println!("map_of_gibbs: {:?}", map_of_gibbs[substance]);
            println!("map_of_gibbs_sym: {:?}", map_of_gibbs_sym[substance]);
        
        }
for substance in &thermo.vec_of_substances {
            let dG_value = map_of_gibbs[substance];

            let dG_from_fun = map_of_gibbs_fun.get(substance).unwrap()(400.0, concentration.clone());
            assert_relative_eq!(dG_value, dG_from_fun, epsilon = 1e-8);

            let dG_sym = map_of_gibbs_sym[substance].clone();
            let dG_from_sym  = dG_sym.lambdify_owned(vec!["T", "w1", "w2"]);
            let dG_from_sym  = dG_from_sym(vec![400.0, 0.5, 0.5]);
            assert_relative_eq!(dG_value, dG_from_sym, epsilon = 1e-8);
        
        }
        

    let _ = thermo.pretty_print();


```

## NIST_scrapper
NIST Chemistry WebBook contains huge amount of thermochemical data. This module can scrap 
data automatically and use it 
```rust
use KiThe::Thermodynamics::DBhandlers::NIST_parser::{ Phase, SearchType, NistParser};
let parser = NistParser::new();

// Example usage
let substance = "CH4";
match parser.get_data(substance, SearchType::All, Phase::Gas) {
    Ok(data) =>{ println!("Data for {}: {:?}", substance, data);
    data.pretty_print();
    #[allow(non_snake_case)]
    let (Cp, dh, ds) = data.caclc_cp_dh_ds(298.15).expect("Error calculating cp, dh, ds");
    println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);},
    Err(e) => eprintln!("Error: {}", e),
                }

let parser = NistParser::new();
let substance = "NaCl";
match parser.get_data(substance, SearchType::All, Phase::Solid) {
    Ok(data) =>{ println!("Data for {}: {:?}", substance, data);
            
    data.pretty_print();
    #[allow(non_snake_case)]
    let (Cp, dh, ds) = data.caclc_cp_dh_ds(298.15).expect("Error calculating cp, dh, ds");
    println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);},
                Err(e) => eprintln!("Error: {}", e),
                
            }
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




