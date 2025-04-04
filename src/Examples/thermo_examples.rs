use crate::Thermodynamics::DBhandlers::CEAdata::CEAdata;
use crate::Thermodynamics::DBhandlers::NASAdata::NASAdata;
use crate::Thermodynamics::DBhandlers::TRANSPORTdata::TransportData;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use RustedSciThe::symbolic::symbolic_engine::Expr;
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
            let _allowed_libs: Vec<String> =
                vec!["Aramco_transport".to_string(), "CEA".to_string()];
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
            let mut tr = TransportData::new();
            let _ = tr.from_serde(CO_data.clone());
            let _ = tr.set_M_unit(Some("g/mol".to_owned()));
            let _ = tr.set_P_unit(Some("atm".to_owned()));
            let _ = tr.set_V_unit(Some("mkPa*s".to_owned()));
            let _ = tr.set_lambda_unit(Some("mW/m/K".to_owned()));
            let T = 473.15; // K 
            tr.P = 1.0;
            tr.M = 28.0; // g/mol
            let _ = tr.calculate_Visc(T);
            assert_relative_eq!(tr.V, 25.2, epsilon = 5.0);
            println!("Viscosity: {:?} mkPa*s", tr.V);

            let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            let mut NASA = NASAdata::new();
            let _ = NASA.from_serde(CO_data.clone());
            let _ = NASA.extract_coefficients(T);
            NASA.calculate_Cp_dH_dS(T);
            let Cp = NASA.Cp;
            println!("Cp: {}", Cp,);
            let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
            let L = tr.calculate_Lambda(Cp, Some(ro), T).unwrap();
            println!("Lambda: {}", L);
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
            // C
            let taylor = NASA.Taylor_series_Cp_dH_dS(300.0, 3);
            let Cp_taylor = taylor.0;
            let dh_taylor = taylor.1;
            let ds_taylor = taylor.2;
            println!(
                "Cp_taylor: {}, Cp exect {}",
                Cp_taylor.lambdify1D()(400.0),
                Cp_value
            );
            println!(
                "dh_taylor: {}, dh exect {}",
                dh_taylor.lambdify1D()(400.0),
                dh_value
            );
            println!(
                "ds_taylor: {}, ds exect {}",
                ds_taylor.lambdify1D()(400.0),
                ds_value
            );
            println!();
        }
        4 => {
            use std::time::Instant;
            let now = Instant::now();
            let thermo_data = ThermoData::new();
            let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            println!("CO_data: {}", CO_data);
            let mut tr = TransportData::new();
            let _ = tr.from_serde(CO_data.clone());
            let _ = tr.set_M_unit(Some("g/mol".to_owned()));
            let _ = tr.set_P_unit(Some("atm".to_owned()));
            let _ = tr.set_V_unit(Some("mkPa*s".to_owned()));
            let _ = tr.set_lambda_unit(Some("mW/m/K".to_owned()));
            let T = 2000.15; // K 
            tr.P = 1.0;
            tr.M = 28.0; // g/mol
            let _ = tr.calculate_Visc(T);

            println!("Viscosity: {:?} mkPa*s", tr.V);

            let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            let mut NASA = NASAdata::new();
            let _ = NASA.from_serde(CO_data.clone());
            let _ = NASA.extract_coefficients(T);
            let _ = NASA.create_sym_Cp_dH_dS();
            let Cp_sym = NASA.clone().Cp_sym;
            NASA.calculate_Cp_dH_dS(T);
            let Cp = NASA.Cp;
            println!("Cp: {}", Cp,);
            let ro = (tr.P * 101325.0) * (tr.M / 1000.0) / (R * T);
            let L = &tr.calculate_Lambda(Cp, Some(ro), T).unwrap();
            let _ = tr.create_closure_Lambda(Cp, Some(ro));
            let Lambda_closure = &mut tr.Lambda_fun;
            let Lambda_from_closure = Lambda_closure(T);

            assert_eq!(Lambda_from_closure, L.clone());
            // let us prove that taylor series of order 2 is enough to calculate Lambda in the temperature range from 300 to 2000
            let taylor_series_Lambda = tr
                .Taylor_series_Lambda(Cp_sym, Expr::Const(ro), 800.0, 2)
                .unwrap();
            let taylor_series_Lambda = taylor_series_Lambda.lambdify1D()(T);
            let elapsed = now.elapsed().as_secs_f64();
            println!("Elapsed: {:.2?}", elapsed);
            println!(
                "Lambda: {:?}, taylor_series_Lambda: {:?}",
                L, taylor_series_Lambda
            );
            assert_relative_eq!(taylor_series_Lambda, L.clone(), epsilon = 3.0);
        }
        5 => {
            use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
            use crate::Thermodynamics::DBhandlers::thermo_api::create_thermal_by_name;
            let thermo_data = ThermoData::new();
            let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
            let CO_data = sublib.get("CO").unwrap();
            println!(" CO data \n {} \n", CO_data);
            let mut NASA = create_thermal_by_name("NASA_gas");
            let _ = NASA.newinstance();
            let _ = NASA.from_serde(CO_data.clone());
            //  assert!(NASA.from_serde(CO_data.clone()).is_ok());
            print!(" this is NASA instance: \n");
            let _ = NASA.print_instance();
            assert!(NASA.extract_model_coefficients(400.0).is_ok());

            let coeffs_len = {
                let coeff_vec = NASA.get_coefficients().unwrap();
                coeff_vec.len()
            };
            assert_eq!(coeffs_len, 7);

            let _ = NASA.calculate_Cp_dH_dS(400.0);
            let Cp = NASA.get_Cp().unwrap();
            let dh = NASA.get_dh().unwrap();
            let ds = NASA.get_ds().unwrap();

            println!("Cp: {}, dh: {}, ds: {}", Cp, dh, ds);
            assert!(Cp > 0.0);
            assert!(dh < 0.0);
            assert!(ds > 0.0);

            let t = 400.0;
            let _ = NASA.create_closures_Cp_dH_dS();

            let Cp_fun = NASA.get_C_fun().unwrap();
            let dh_fun = NASA.get_dh_fun().unwrap();
            let ds_fun = NASA.get_ds_fun().unwrap();
            assert_relative_eq!((Cp_fun)(t), Cp, epsilon = 1e-6);
            assert_relative_eq!((dh_fun)(t), dh, epsilon = 1e-6);
            assert_relative_eq!((ds_fun)(t), ds, epsilon = 1e-6);

            let _ = NASA.create_sym_Cp_dH_dS();
            let Cp_sym = NASA.get_Cp_sym().unwrap();
            let Cp_T = Cp_sym.lambdify1D();
            let Cp_value = Cp_T(400.0);
            assert_relative_eq!(Cp_value, Cp, epsilon = 1e-6);

            let dh_sym = NASA.get_dh_sym().unwrap();
            let dh_T = dh_sym.lambdify1D();
            let dh_value = dh_T(400.0);
            assert_relative_eq!(dh_value, dh, epsilon = 1e-6);

            let ds_sym = NASA.get_ds_sym().unwrap();
            let ds_T = ds_sym.lambdify1D();
            let ds_value = ds_T(400.0);
            assert_relative_eq!(ds_value, ds, epsilon = 1e-6);
            _ = NASA.pretty_print_data();
        }
        _ => println!("Invalid task number"),
    }
}
