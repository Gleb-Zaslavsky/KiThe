/// module for classical thermodynamics and chemical equilibrium
pub mod ChemEquilibrium;
#[allow(non_snake_case)]
/// handlers for different formats of thermodynamics and heat-mass transfer data
pub mod DBhandlers;
/// heat-mass transfer data agregator
pub mod User_substances;
/// tests
pub mod User_substances_tests;
/// main functionality to open thermodynamics and heat-mass transfer libraries
pub mod thermo_lib_api;
/// calculations of thermodynamic properties, creating closures and symbolic expressions
pub mod thermo_properties_api;
/// calculations of heat-mass transfer properties, creating closures and symbolic expressions
pub mod transport_properties_api;
