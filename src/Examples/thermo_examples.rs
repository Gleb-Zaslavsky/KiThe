use crate::Thermodynamics::DBhandlers::CEAdata::CEAdata;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::Thermodynamics::DBhandlers::NASAdata::NASAdata;
use crate::Thermodynamics::DBhandlers::TRANSPORTdata::TransportData;
use approx::assert_relative_eq;
const R: f64 = 8.314;
pub fn thermo_examples(thermotask: usize) {
    //

    match thermotask {
        0 => {
            let mut thermo_data = ThermoData::new();
            println!("Libraries on board {:?} \n \n", thermo_data.AllLibraries);
            let this_lib_subs = thermo_data.subs_of_this_lib("Cantera_nasa_base_gas");
            println!("subs of this lib {:?} \n \n", this_lib_subs);
            let allowed_libs: Vec<String> = vec!["Aramco_transport".to_string(), "CEA".to_string()];
            thermo_data.search_libs_for_subs(vec!["CO".to_string()], None);
            println!(
                "hashmap_of_thermo_data {:?}",
                thermo_data.hashmap_of_thermo_data
            );
            thermo_data.pretty_print_thermo_data();
        }
     
        1 => {

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
        }
        2 => {
            let thermo_data = ThermoData::new();
            let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            println!("CO_data: {}", CO_data);
            let mut tr = TransportData ::new();
            tr.from_serde(CO_data.clone());
            tr.set_M_unit(Some("g/mol".to_owned())  );
            tr.set_P_unit(Some("atm".to_owned()));
            tr.set_V_unit( Some("mkPa*s".to_owned()));
            tr. set_lambda_unit(Some("mW/m/K".to_owned()));
            let T = 473.15; // K 
            tr.P = 1.0;
            tr.M = 28.0; // g/mol
            tr.calculate_Visc(T);
            assert_relative_eq!( tr.V, 25.2, epsilon = 5.0);
            println!("Viscosity: {:?} mkPa*s", tr.V);
            let P = 1.0; // atm
        
            let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            let mut NASA = NASAdata::new();
            NASA.from_serde(CO_data.clone());
            NASA.extract_coefficients(T);
            NASA.calculate_Cp_dH_dS(T);
            let Cp = NASA.Cp;
            println!("Cp: {}", Cp, );
            let ro = (tr.P*101325.0)*(tr.M/1000.0)/(R*T);
            let L = tr.calculate_Lambda(Cp, ro, T);
            println!("Lambda: {}",L);
    
        }
        3 => {
            let thermo_data = ThermoData::new();
            let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            let mut NASA = NASAdata::new();
    
            assert!(NASA.from_serde(CO_data.clone()).is_ok());
            assert!(NASA.extract_coefficients(400.0).is_ok());
    
            println!("\n \n {:?} \n \n ", NASA.coeffs);
            let coeffs_len = {
                let (_a, _b, _c, _d, _e, _f, _g) = NASA.coeffs;
                7 // Since we know the tuple has 7 elements
            };
            assert_eq!(coeffs_len, 7);
    
            NASA.calculate_Cp_dH_dS(400.0);
            let Cp = NASA.Cp;
            let dh = NASA.dh;
            let ds = NASA.ds;
    
            println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
            assert!(Cp > 0.0);
            assert!(dh < 0.0);
            assert!(ds > 0.0);
    
            let t = 400.0;
            NASA.create_closures_Cp_dH_dS();
    
            let Cp_fun = &NASA.C_fun;
            let dh_fun = &NASA.dh_fun;
            let ds_fun = &NASA.ds_fun;
            assert_relative_eq!((Cp_fun)(t), NASA.Cp, epsilon = 1e-6);
            assert_relative_eq!((dh_fun)(t), NASA.dh, epsilon = 1e-6);
            assert_relative_eq!((ds_fun)(t), NASA.ds, epsilon = 1e-6);
    
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
        }
        _ => println!("Invalid task number"),
    }
}
