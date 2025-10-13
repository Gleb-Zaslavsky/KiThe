///handler of NASA format of thermodynamic data
/// # Examples
///
/// ```
/// use KiThe::Thermodynamics::thermo_lib_api::ThermoData;
/// use KiThe::Thermodynamics::DBhandlers::NASAdata::NASAdata;
/// use approx::assert_relative_eq;
///  // open a database
/// let thermo_data = ThermoData::new();
/// // we a interested in NASA_gas library
/// let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
/// // we are interested in CO substance
/// let CO_data = sublib.get("CO").unwrap();
/// // we create a new NASAdata object
/// let mut NASA = NASAdata::new();
/// // we fill it with data
/// NASA.from_serde(CO_data.clone());
/// // get the coefficients from the certtain temperature (400 K)
/// NASA.extract_coefficients(400.0);
/// println!("\n \n {:?} \n \n ", NASA.coeffs);
/// let coeffs_len = {
///    let (_a, _b, _c, _d, _e, _f, _g) = NASA.coeffs;
///    7 // Since we know the tuple has 7 elements
/// };
/// assert_eq!(coeffs_len, 7);
/// // calculate values Cp, dH, dS at 400 K
/// NASA.calculate_Cp_dH_dS(400.0);
/// let Cp = NASA.Cp;
/// let dh = NASA.dh;
/// let ds = NASA.ds;
/// println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
/// let t = 400.0;
/// // create closures Box<dyn Fn(f64)- f64> for Cp, dH, dS
/// NASA.create_closures_Cp_dH_dS();
/// let Cp_fun = &NASA.C_fun;
/// let dh_fun = &NASA.dh_fun;
/// let ds_fun = &NASA.ds_fun;
/// assert_relative_eq!((Cp_fun)(t), NASA.Cp, epsilon = 1e-6);
/// assert_relative_eq!((dh_fun)(t), NASA.dh, epsilon = 1e-6);
/// assert_relative_eq!((ds_fun)(t), NASA.ds, epsilon = 1e-6);
/// // create symbolic expressions for Cp, dH, dS
/// NASA.create_sym_Cp_dH_dS();
/// let Cp_sym = &NASA.Cp_sym;
/// // turn symbolic expressions into functions
/// let Cp_T = Cp_sym.lambdify1D();
///let Cp_value = Cp_T(400.0);
/// assert_relative_eq!(Cp_value, NASA.Cp, epsilon = 1e-6);
/// let dh_sym = &NASA.dh_sym;
/// let dh_T = dh_sym.lambdify1D();
/// let dh_value = dh_T(400.0);
/// assert_relative_eq!(dh_value, NASA.dh, epsilon = 1e-6);
/// let ds_sym = &NASA.ds_sym;
/// let ds_T = ds_sym.lambdify1D();
/// let ds_value = ds_T(400.0);
/// assert_relative_eq!(ds_value, NASA.ds, epsilon = 1e-6);
///  // calculate Taylor series for Cp, dH, dS at 300 K and 3rd order
/// let taylor = NASA.Taylor_series_Cp_dH_dS(300.0, 3);
/// let Cp_taylor = taylor.0;
/// let dh_taylor = taylor.1;
/// let ds_taylor = taylor.2;
/// println!(
///    "Cp_taylor: {}, Cp exect {}",
///    Cp_taylor.lambdify1D()(400.0),
///    Cp_value
/// );
/// println!(
///    "dh_taylor: {}, dh exect {}",
///    dh_taylor.lambdify1D()(400.0),
///    dh_value
///);
/// println!(
///    "ds_taylor: {}, ds exect {}",
///    ds_taylor.lambdify1D()(400.0),
///    ds_value
///);
/// ```
pub mod NASAdata;

/// heat and mass transfer data from CEA (Chemical Equilibrium with Applications)
/// # Examples
/// ```
/// use KiThe::Thermodynamics::thermo_lib_api::ThermoData;
/// use KiThe::Thermodynamics::DBhandlers::CEAdata::CEAdata;
/// use approx::assert_relative_eq;
/// // open a database
/// let thermo_data = ThermoData::new();
/// // we a interested in CEA library
/// let sublib = thermo_data.LibThermoData.get("CEA").unwrap();
/// // we are interested in CO substance
/// let CO_data = sublib.get("CO").unwrap();
/// // we create a new CEAdata object
/// let mut CEA = CEAdata::new();
/// // we fill it with data
/// CEA.from_serde(CO_data.clone()).unwrap();
/// // set units
/// CEA.set_lambda_unit(Some("mW/m/K".to_string())).unwrap();
/// CEA.set_V_unit(Some("mkPa*s".to_string())).unwrap();
/// // parse coefficients for Lambda and Viscosity at 500 K
/// CEA.parse_coefficients().unwrap();
/// // calculate Lambda and Viscosity at 500 K
/// CEA.extract_coefficients(500.0).unwrap();
/// // calculate Lambda and Viscosity at 500 K
/// let lambda = CEA.calculate_Lambda(500.0).unwrap();
/// let visc = CEA.calculate_Visc(500.0).unwrap();
/// println!("Lambda, mW/m/K: {:?}, Visc: {:?}", lambda, visc);
/// // create closures for Lambda and Viscosity
/// let lambda_closure = CEA.create_closure_Lambda().unwrap();
///  let visc_closure = CEA.create_closure_Visc().unwrap();
/// // calculate Lambda and Viscosity at 500 K
/// let lamda_value = lambda_closure (500.0);
/// let visc_value = visc_closure(500.0);
/// // compare with directly calculated values
/// assert_relative_eq!(lambda, lamda_value, epsilon = 1e-6);
/// assert_relative_eq!(visc, visc_value, epsilon = 1e-6);
/// // create symbolic expressions for Lambda and Viscosity
/// assert!(CEA.create_sym_Lambda().is_ok());
/// assert!(CEA.create_sym_Visc().is_ok());
/// // turn symbolic expressions into functions
/// let lambda_sym = CEA.Lambda_sym.as_ref().unwrap();
/// let lambda_val = lambda_sym.lambdify1D()(500.0);
/// let visc_sym = CEA.V_sym.as_ref().unwrap();
/// let visc_val = visc_sym.lambdify1D()(500.0);
/// // compare with directly calculated values
/// assert_relative_eq!(lambda_val, lambda, epsilon = 1e-6);
/// assert_relative_eq!(visc_val, visc, epsilon = 1e-6);
/// ```
pub mod CEAdata;
/// parser of NIST Chemical WebBook Cp, dH, dS data
/// # Examples
/// ```
/// use KiThe::Thermodynamics::DBhandlers::NIST_parser::NistParser;
/// use KiThe::Thermodynamics::DBhandlers::NIST_parser::SearchType;
/// use KiThe::Thermodynamics::DBhandlers::NIST_parser::Phase;
///
/// let parser = NistParser::new();
/// // Example usage
/// let parser = NistParser::new();
/// let substance = "NaCl";
/// match parser.get_data(substance, SearchType::All, Phase::Solid) {
///    Ok(mut data) => {
///        println!("Data for {}: {:?}", substance, data);
///        data.pretty_print();
///         let _ = data.extract_coefficients(298.15);
///        #[allow(non_snake_case)]
///        let (Cp, dh, ds) = data
///            .caclc_cp_dh_ds(298.15)
///            .expect("Error calculating cp, dh, ds");
///        println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);
///    }
///    Err(e) => eprintln!("Error: {}", e),
///}
/// match parser.get_data(substance, SearchType::All, Phase::Liquid) {
///    Ok(mut data) => {
///        println!("Data for {}: {:?}", substance, data);
///        data.pretty_print();
///        let _ = data.extract_coefficients(1200.15);
///        #[allow(non_snake_case)]
///        let (Cp, dh, ds) = data
///            .caclc_cp_dh_ds(1200.15)
///            .expect("Error calculating cp, dh, ds");
///        println!("Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}", Cp, dh, ds);
///    }
///    Err(e) => eprintln!("Error: {}", e),
/// }
///    
/// ```
pub mod NIST_parser;
/// handler of NIST format of thermodynamic data
/// # Examples
/// ```
/// use KiThe::Thermodynamics::DBhandlers::NIST_parser::SearchType;
/// use KiThe::Thermodynamics::DBhandlers::NIST_parser::Phase;
/// use KiThe::Thermodynamics::DBhandlers::NISTdata::NISTdata;
/// use approx::assert_relative_eq;
///         let mut nist = NISTdata::new();
/// // scrape NIST site for CO data
/// let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);
/// // Test calculate_Cp_dH_dS
/// let _ = nist.extract_coefficients(400.0);
/// nist.calculate_Cp_dH_dS(400.0);
/// assert!(nist.Cp > 0.0);
/// assert!(nist.dh != 0.0);
/// assert!(nist.ds != 0.0);

/// // Test create_closures_Cp_dH_dS
/// nist.create_closures_Cp_dH_dS();
/// let t = 400.0;
/// assert_relative_eq!((nist.C_fun)(t), nist.Cp, epsilon = 1e-6);
/// assert_relative_eq!((nist.dh_fun)(t), nist.dh, epsilon = 1e-6);
/// assert_relative_eq!((nist.ds_fun)(t), nist.ds, epsilon = 1e-6);
///
/// // Test create_sym_Cp_dH_dS
/// nist.create_sym_Cp_dH_dS();
/// let Cp_sym = &nist.Cp_sym;
/// let Cp_T = Cp_sym.lambdify1D();
/// let Cp_value = Cp_T(400.0);
/// assert_relative_eq!(Cp_value, nist.Cp, epsilon = 1e-6);
// Test pretty_print_data
/// assert!(nist.pretty_print().is_ok());
/// /// // Test Taylor_series_cp_dh_ds
/// let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_Cp_dH_dS(400.0, 3).unwrap();
/// assert!(!Cp_taylor.is_zero());
/// assert!(!dh_taylor.is_zero());
/// assert!(!ds_taylor.is_zero());
///
/// ```
pub mod NISTdata;
/// heat and mass transfer data based on collision integrals
pub mod TRANSPORTdata;
/// general api for thermodynamic data (NASA, NIST)
///
/// general api for transport (heat and mass transfer) data (CEA, TRANSPORT)
/// # Examples
/// ```
/// use KiThe::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
/// use KiThe::Thermodynamics::DBhandlers::thermo_api::create_thermal_by_name;
/// use  KiThe::Thermodynamics::thermo_lib_api::ThermoData;
/// use approx::assert_relative_eq;
/// let thermo_data = ThermoData::new();
/// let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
/// let CO_data = sublib.get("CO").unwrap();
/// println!(" CO data \n {} \n", CO_data);
/// let mut NASA = create_thermal_by_name("NASA_gas");
/// let _ = NASA.newinstance();
/// let _ = NASA.from_serde(CO_data.clone());
/// //  assert!(NASA.from_serde(CO_data.clone()).is_ok());
/// print!(" this is NASA instance: \n");
/// let _ = NASA.print_instance();
/// assert!(NASA.extract_model_coefficients(400.0).is_ok());
///
/// let coeffs_len = {
///    let coeff_vec = NASA.get_coefficients().unwrap();
///    coeff_vec.len()
/// };
/// assert_eq!(coeffs_len, 7);
///
/// let _ = NASA.calculate_Cp_dH_dS(400.0);
/// let Cp = NASA.get_Cp().unwrap();
/// let dh = NASA.get_dh().unwrap();
/// let ds = NASA.get_ds().unwrap();
///
/// println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
/// assert!(Cp > 0.0);
/// assert!(dh < 0.0);
/// assert!(ds > 0.0);
///
/// let t = 400.0;
/// let _ = NASA.create_closures_Cp_dH_dS();
///
/// let Cp_fun = NASA.get_C_fun().unwrap();
/// let dh_fun = NASA.get_dh_fun().unwrap();
/// let ds_fun = NASA.get_ds_fun().unwrap();
/// assert_relative_eq!((Cp_fun)(t), Cp, epsilon = 1e-6);
/// assert_relative_eq!((dh_fun)(t), dh, epsilon = 1e-6);
/// assert_relative_eq!((ds_fun)(t), ds, epsilon = 1e-6);
///
/// let _ = NASA.create_sym_Cp_dH_dS();
/// let Cp_sym = NASA.get_Cp_sym().unwrap();
/// let Cp_T = Cp_sym.lambdify1D();
/// let Cp_value = Cp_T(400.0);
/// assert_relative_eq!(Cp_value, Cp, epsilon = 1e-6);
///
/// let dh_sym = NASA.get_dh_sym().unwrap();
/// let dh_T = dh_sym.lambdify1D();
/// let dh_value = dh_T(400.0);
/// assert_relative_eq!(dh_value, dh, epsilon = 1e-6);
///
/// let ds_sym = NASA.get_ds_sym().unwrap();
/// let ds_T = ds_sym.lambdify1D();
/// let ds_value = ds_T(400.0);
/// assert_relative_eq!(ds_value, ds, epsilon = 1e-6);
/// _ = NASA.pretty_print_data();
/// ```
pub mod thermo_api;

pub mod transport_api;
//pub use transport_api::{TransportCalculator, TransportError, LambdaUnit, ViscosityUnit};
mod NASAdata_tests;
mod NIST_parser_tests;
