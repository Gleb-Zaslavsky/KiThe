///////////////////////////TESTING////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use std::result;

    use crate::Thermodynamics::DBhandlers::NASAdata::{Coeffs, Cp, NASAError, NASAdata, dh, ds};
    use crate::Thermodynamics::DBhandlers::thermo_api::{
        EnergyUnit, ThermoCalculator, energy_dimension,
    };
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use approx::assert_relative_eq;

    #[test]
    fn test_new() {
        let nasa_data = NASAdata::new();
        assert_eq!(nasa_data.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
        assert_eq!(nasa_data.Cp, 0.0);
        assert_eq!(nasa_data.dh, 0.0);
        assert_eq!(nasa_data.ds, 0.0);
    }

    #[test]
    fn test_extract_coefficients() {
        let mut nasa_data = NASAdata::new();
        nasa_data.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        assert!(nasa_data.extract_coefficients(500.0).is_ok());
        assert_eq!(nasa_data.coeffs, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0));

        assert!(nasa_data.extract_coefficients(700.0).is_ok());
        assert_eq!(nasa_data.coeffs, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0));
    }

    #[test]
    fn test_extract_coefficients_out_of_range() {
        let mut nasa_data = NASAdata::new();
        nasa_data.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];

        let result = nasa_data.extract_coefficients(1500.0);
        assert!(matches!(
            result,
            Err(NASAError::NoCoefficientsFound {
                temperature: 1500.0,
                ..
            })
        ));
    }

    #[test]
    fn test_set_unit() {
        let mut nasa_data = NASAdata::new();

        assert!(nasa_data.set_unit("J").is_ok());
        assert_eq!(nasa_data.unit, Some("J".to_string()));
        assert_eq!(nasa_data.unit_multiplier, 4.184);

        assert!(nasa_data.set_unit("cal").is_ok());
        assert_eq!(nasa_data.unit, Some("cal".to_string()));
        assert_eq!(nasa_data.unit_multiplier, 1.0);

        let result = nasa_data.set_unit("invalid");
        assert!(matches!(result, Err(NASAError::UnsupportedUnit(unit)) if unit == "invalid"));
    }

    #[test]
    fn test_create_closures_cp_dh_ds() {
        let mut nasa_data = NASAdata::new();
        nasa_data.coeffs = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0);
        let _ = nasa_data.set_unit("cal");
        nasa_data.create_closures_Cp_dH_dS();

        let t = 500.0;
        assert_relative_eq!(
            (nasa_data.C_fun)(t),
            Cp(t, 1.0, 2.0, 3.0, 4.0, 5.0),
            epsilon = 1e-6
        );
        assert_relative_eq!(
            (nasa_data.dh_fun)(t),
            dh(t, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
            epsilon = 1e-6
        );
        assert_relative_eq!(
            (nasa_data.ds_fun)(t),
            ds(t, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0),
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_with_real_data() {
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

    #[test]
    fn test_thermo_calculator_error_handling() {
        let mut nasa = NASAdata::new();

        // Test invalid temperature range
        let result = nasa.extract_model_coefficients(1000.0);
        assert!(result.is_err());

        // Test invalid unit
        let result = nasa.set_unit(&energy_dimension(EnergyUnit::J));
        assert!(result.is_ok());

        // Test invalid serde data
        let invalid_data = serde_json::json!({
            "invalid": "data"
        });
        let result = nasa.from_serde(invalid_data);
        assert!(result.is_err());
    }
    #[test]
    fn test_with_real_data_ThermoCalc() {
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
    }

    #[test]
    fn test_thermo_calculator_nasa() {
        // use super::EnergyUnit;
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());

        // Test set_unit
        // assert!(nasa.set_unit(Some(EnergyUnit::J)).is_ok());
        // assert!(nasa.set_unit(Some(EnergyUnit::Cal)).is_ok());

        // Test extract_model_coefficients
        assert!(nasa.extract_model_coefficients(400.0).is_ok());

        // Test calculate_Cp_dH_dS
        nasa.calculate_Cp_dH_dS(400.0);
        assert!(nasa.Cp > 0.0);
        assert!(nasa.dh != 0.0);
        assert!(nasa.ds != 0.0);

        // Test create_closures_Cp_dH_dS
        nasa.create_closures_Cp_dH_dS();
        let t = 400.0;
        assert_relative_eq!((nasa.C_fun)(t), nasa.Cp, epsilon = 1e-6);
        assert_relative_eq!((nasa.dh_fun)(t), nasa.dh, epsilon = 1e-6);
        assert_relative_eq!((nasa.ds_fun)(t), nasa.ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        nasa.create_sym_Cp_dH_dS();
        let Cp_sym = &nasa.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        let dh_sym = &nasa.dh_sym;
        let dh_T = dh_sym.lambdify1D();
        let dh_value = dh_T(400.0);
        let ds_sym = &nasa.ds_sym;
        let ds_T = ds_sym.lambdify1D();
        let ds_value = ds_T(400.0);
        assert_relative_eq!(Cp_value, nasa.Cp, epsilon = 1e-6);
        assert_relative_eq!(dh_value, nasa.dh, epsilon = 1e-6);
        assert_relative_eq!(ds_value, nasa.ds, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nasa.Taylor_series_cp_dh_ds(300.0, 3).unwrap();
        let Cp_taylor_T = Cp_taylor.lambdify1D();
        let Cp_value = Cp_taylor_T(400.0);
        assert_relative_eq!(Cp_value, nasa.Cp, epsilon = 3.0);
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());

        // Test pretty_print_data
        assert!(nasa.pretty_print_data().is_ok());
    }
    //////////////////////DUMMY////////////////////////////////////////////
    #[test]
    fn test_fitting_adjacent_with_dummy() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(800.0).unwrap();
        assert!(*i == 0);
        let (i, _) = &nasa.interval_for_this_T(1200.0).unwrap();
        assert!(*i == 1);

        // making dummy structsa with identical coefficients
        let coeffs_min_dummy = Coeffs {
            T: (300.0, 500.0),
            coeff: coeffs_min.coeff.clone(),
        };
        let coeffs_max_dummy = Coeffs {
            T: (500.0, 1000.0),
            coeff: coeffs_min.coeff.clone(),
        };
        //  assert!(nasa.T_interval.is_some());
        let result = nasa.fitting_adjacent2(
            coeffs_min_dummy.clone(),
            coeffs_max_dummy.clone(),
            500.0,
            800.0,
        );
        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        assert_relative_eq!(coeffs_result.0, coeffs_min.coeff.0, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.1, coeffs_min.coeff.1, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.2, coeffs_min.coeff.2, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.3, coeffs_min.coeff.3, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.4, coeffs_min.coeff.4, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.5, coeffs_min.coeff.5, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.6, coeffs_min.coeff.6, epsilon = 1e-3);
    }

    #[test]
    fn test_fitting_adjacent_linear_method_with_dummy() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(800.0).unwrap();
        assert!(*i == 0);
        let (i, _) = &nasa.interval_for_this_T(1200.0).unwrap();
        assert!(*i == 1);

        // making dummy structsa with identical coefficients
        let coeffs_min_dummy = Coeffs {
            T: (300.0, 500.0),
            coeff: coeffs_min.coeff.clone(),
        };
        let coeffs_max_dummy = Coeffs {
            T: (500.0, 1000.0),
            coeff: coeffs_min.coeff.clone(),
        };
        //  assert!(nasa.T_interval.is_some());
        let result = nasa.fit_cp_dh_ds(
            coeffs_min_dummy.clone(),
            coeffs_max_dummy.clone(),
            500.0,
            800.0,
        );
        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        assert_relative_eq!(coeffs_result.0, coeffs_min.coeff.0, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.1, coeffs_min.coeff.1, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.2, coeffs_min.coeff.2, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.3, coeffs_min.coeff.3, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.4, coeffs_min.coeff.4, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.5, coeffs_min.coeff.5, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.6, coeffs_min.coeff.6, epsilon = 1e-3);
    }
    #[test]
    fn test_fitting_adjacent3_with_dummy() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(800.0).unwrap();
        assert!(*i == 0);
        let (i, _) = &nasa.interval_for_this_T(1200.0).unwrap();
        assert!(*i == 1);

        // making dummy structsa with identical coefficients
        let coeffs_min_dummy = Coeffs {
            T: (300.0, 500.0),
            coeff: coeffs_min.coeff.clone(),
        };
        let coeffs_max_dummy = Coeffs {
            T: (500.0, 1000.0),
            coeff: coeffs_min.coeff.clone(),
        };
        //  assert!(nasa.T_interval.is_some());
        let result = nasa.fitting_adjacent3(
            coeffs_min_dummy.clone(),
            coeffs_max_dummy.clone(),
            500.0,
            800.0,
        );
        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        assert_relative_eq!(coeffs_result.0, coeffs_min.coeff.0, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.1, coeffs_min.coeff.1, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.2, coeffs_min.coeff.2, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.3, coeffs_min.coeff.3, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.4, coeffs_min.coeff.4, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.5, coeffs_min.coeff.5, epsilon = 1e-3);
        assert_relative_eq!(coeffs_result.6, coeffs_min.coeff.6, epsilon = 1e-3);
    }
    //////////////////////REAL DATA/////////////////////////////////////////////////////
    #[test]
    fn test_fitting_adjacent_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        nasa.set_T_interval(900.0, 1100.0);
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(900.0).unwrap();
        assert!(*i == 0);
        let (i, coeffs_max) = &nasa.interval_for_this_T(1100.0).unwrap();
        assert!(*i == 1);
        let g_min = coeffs_min.coeff.6;
        let g_max = coeffs_max.coeff.6;
        println!("g_min {:?} \n", g_min);
        println!("g_max {:?} \n", g_max);
        assert!(nasa.T_interval.is_some());
        let result = nasa.fitting_adjacent2(coeffs_min.clone(), coeffs_max.clone(), 900.0, 1100.0);

        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        println!("coeffs for  min T {:?} \n", coeffs_min.coeff);
        println!("coeffs for  max T {:?} \n", coeffs_max.coeff);
        println!("{:?}", result.unwrap());
    }

    #[test]
    fn test_fitting_adjacent3_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        nasa.set_T_interval(900.0, 1100.0);
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(900.0).unwrap();
        assert!(*i == 0);
        let (i, coeffs_max) = &nasa.interval_for_this_T(1100.0).unwrap();
        assert!(*i == 1);
        let g_min = coeffs_min.coeff.6;
        let g_max = coeffs_max.coeff.6;
        println!("g_min {:?} \n", g_min);
        println!("g_max {:?} \n", g_max);
        assert!(nasa.T_interval.is_some());
        let result = nasa.fitting_adjacent3(coeffs_min.clone(), coeffs_max.clone(), 900.0, 1100.0);

        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        println!("coeffs for  min T {:?} \n", coeffs_min.coeff);
        println!("coeffs for  max T {:?} \n", coeffs_max.coeff);
        println!("{:?}", result.unwrap());
    }

    #[test]
    fn test_fitting_adjacent_linear_method_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        // Test newinstance
        assert!(nasa.newinstance().is_ok());

        // Test from_serde
        assert!(nasa.from_serde(CO_data.clone()).is_ok());
        nasa.set_T_interval(900.0, 1100.0);
        println!(" this is NASA instance: {:?} \n", nasa);
        nasa.parse_coefficients().unwrap();
        assert!(!nasa.coeffs_map.is_empty());
        println!("coeffs {:?} \n", nasa.coeffs_map);
        let (i, coeffs_min) = &nasa.interval_for_this_T(900.0).unwrap();
        assert!(*i == 0);
        let (i, coeffs_max) = &nasa.interval_for_this_T(1100.0).unwrap();
        assert!(*i == 1);
        let g_min = coeffs_min.coeff.6;
        let g_max = coeffs_max.coeff.6;
        println!("g_min {:?} \n", g_min);
        println!("g_max {:?} \n", g_max);
        assert!(nasa.T_interval.is_some());
        let result = nasa.fit_cp_dh_ds(coeffs_min.clone(), coeffs_max.clone(), 900.0, 1100.0);

        assert!(result.is_ok());
        let coeffs_result = nasa.coeffs;
        println!("coeffs result {:?} \n", coeffs_result);
        println!("coeffs for  min T {:?} \n", coeffs_min.coeff);
        println!("coeffs for  max T {:?} \n", coeffs_max.coeff);
        println!("{:?}", result.unwrap());
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////
    #[test]
    fn test_fitting_adjacent_real_data_for_several_subs() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let subs = vec!["CO", "CO2", "NO", "N2"];
        for sub in subs {
            let CO_data = sublib.get(sub).unwrap();
            let mut nasa = NASAdata::new();

            // Test newinstance
            assert!(nasa.newinstance().is_ok());

            // Test from_serde
            assert!(nasa.from_serde(CO_data.clone()).is_ok());
            nasa.set_T_interval(900.0, 1100.0);
            // println!(" this is NASA instance: {:?} \n", nasa);
            nasa.parse_coefficients().unwrap();
            assert!(!nasa.coeffs_map.is_empty());

            let data_for_interva0 = nasa.coeffs_map.get(&(0 as usize)).unwrap();
            // creating interval
            let T_ceter1 = data_for_interva0.T.1;
            let data_for_interva1 = nasa.coeffs_map.get(&(1 as usize)).unwrap();
            let T_ceter2 = data_for_interva1.T.0;
            assert_eq!(T_ceter1, T_ceter2);
            let T_min = T_ceter1 - 500.0;
            let T_max = T_ceter1 + 500.0;
            nasa.set_T_interval(T_min, T_max);
            let res = nasa.fitting_coeffs_for_T_interval();

            assert!(res.is_ok());
        }
    }

    #[test]
    fn test_fitting_adjacent_real_data_for_several_subs2() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let subs = vec!["NH3", "O2", "N2"];
        for sub in subs {
            let CO_data = sublib.get(sub).unwrap();
            let mut nasa = NASAdata::new();

            // Test newinstance
            assert!(nasa.newinstance().is_ok());

            // Test from_serde
            assert!(nasa.from_serde(CO_data.clone()).is_ok());
            nasa.set_T_interval(900.0, 1100.0);
            // println!(" this is NASA instance: {:?} \n", nasa);
            nasa.parse_coefficients().unwrap();
            assert!(!nasa.coeffs_map.is_empty());

            let data_for_interva0 = nasa.coeffs_map.get(&(0 as usize)).unwrap();
            // creating interval
            let T_ceter1 = data_for_interva0.T.1;
            let data_for_interva1 = nasa.coeffs_map.get(&(1 as usize)).unwrap();
            let T_ceter2 = data_for_interva1.T.0;
            assert_eq!(T_ceter1, T_ceter2);
            let T_min = T_ceter1 - 500.0;
            let T_max = T_ceter1 + 500.0;
            nasa.set_T_interval(T_min, T_max);
            let res = nasa.fitting_coeffs_for_T_interval();

            assert!(res.is_ok());
        }
    }
    #[test]
    fn test_fitting_coeffs_multiple_substances() {
        let thermo_data = ThermoData::new();
        let subs = vec!["CO".to_string(), "CO2".to_string(), "N2".to_string()];
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        for sub in subs {
            let CO_data = sublib.get(&sub).unwrap();
            let mut nasa = NASAdata::new();

            // Test newinstance
            assert!(nasa.newinstance().is_ok());

            // Test from_serde
            assert!(nasa.from_serde(CO_data.clone()).is_ok());
            nasa.set_T_interval(600.0, 1000.0);
            // println!(" this is NASA instance: {:?} \n", nasa);
            nasa.parse_coefficients().unwrap();
            assert!(!nasa.coeffs_map.is_empty());

            let res = nasa.fitting_coeffs_for_T_interval();

            assert!(res.is_ok());
        }
    }

    #[test]
    fn test_interval_for_this_T() {
        let mut nasa = NASAdata::new();
        nasa.input.Cp = vec![
            300.0, 1000.0, 3000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
            8.0,
        ];
        nasa.parse_coefficients().unwrap();

        let (interval_idx, coeffs) = nasa.interval_for_this_T(500.0).unwrap();
        assert_eq!(interval_idx, 0);
        assert_eq!(coeffs.T, (300.0, 1000.0));

        let (interval_idx, coeffs) = nasa.interval_for_this_T(1500.0).unwrap();
        assert_eq!(interval_idx, 1);
        assert_eq!(coeffs.T, (1000.0, 3000.0));
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_same_interval() {
        let mut nasa = NASAdata::new();
        nasa.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        nasa.parse_coefficients().unwrap();
        nasa.set_T_interval(400.0, 800.0); // Both in same interval

        assert!(nasa.fitting_coeffs_for_T_interval().is_ok());
        assert_eq!(nasa.coeffs, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0));
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_no_interval_set() {
        let mut nasa = NASAdata::new();

        assert!(matches!(
            nasa.fitting_coeffs_for_T_interval(),
            Err(NASAError::InvalidTemperatureRange)
        ));
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_adjacent_intervals() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();
        nasa.set_T_interval(900.0, 1100.0); // Adjacent intervals
        nasa.norm_threshold = 1.0; // Relax threshold for test

        assert!(nasa.fitting_coeffs_for_T_interval().is_ok());
        // Coeffs should be fitted values, not original
        assert_ne!(nasa.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_non_adjacent() {
        let mut nasa = NASAdata::new();
        nasa.input.Cp = vec![
            300.0, 1000.0, 2000.0, 3000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0, 4.0, 5.0,
            6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
        ];
        nasa.parse_coefficients().unwrap();
        nasa.set_T_interval(400.0, 1500.0); // Non-adjacent intervals
        let result = nasa.fitting_coeffs_for_T_interval();
        println!("Result: {:?}", result);
    }

    #[test]
    fn test_fitted_coeffs_closures_and_symbolic() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();
        nasa.set_T_interval(900.0, 1100.0);
        nasa.norm_threshold = 1.0;

        // Fit coefficients for temperature interval
        assert!(nasa.fitting_coeffs_for_T_interval().is_ok());
        let fitted_coeffs = nasa.coeffs;

        // Create closures using fitted coefficients
        nasa.create_closures_Cp_dH_dS();

        // Create symbolic expressions using fitted coeffs
        nasa.create_sym_Cp_dH_dS();

        let test_temp = 1000.0;

        // Calculate values using fitted coefficients
        nasa.calculate_Cp_dH_dS(test_temp);
        let cp_calc = nasa.Cp;
        let dh_calc = nasa.dh;
        let ds_calc = nasa.ds;

        // Test closures use fitted coefficients
        let cp_closure = (nasa.C_fun)(test_temp);
        let dh_closure = (nasa.dh_fun)(test_temp);
        let ds_closure = (nasa.ds_fun)(test_temp);

        assert_relative_eq!(cp_closure, cp_calc, epsilon = 1e-6);
        assert_relative_eq!(dh_closure, dh_calc, epsilon = 1e-6);
        assert_relative_eq!(ds_closure, ds_calc, epsilon = 1e-6);

        // Test symbolic expressions use fitted coefficients
        let cp_sym_func = nasa.Cp_sym.lambdify1D();
        let dh_sym_func = nasa.dh_sym.lambdify1D();
        let ds_sym_func = nasa.ds_sym.lambdify1D();

        let cp_sym_val = cp_sym_func(test_temp);
        let dh_sym_val = dh_sym_func(test_temp);
        let ds_sym_val = ds_sym_func(test_temp);

        assert_relative_eq!(cp_sym_val, cp_calc, epsilon = 1e-6);
        assert_relative_eq!(dh_sym_val, dh_calc, epsilon = 1e-6);
        assert_relative_eq!(ds_sym_val, ds_calc, epsilon = 1e-6);

        // Verify fitted coefficients are different from original
        let original_coeffs = nasa.coeffs_map.get(&0).unwrap().coeff;
        assert_ne!(fitted_coeffs, original_coeffs);
    }

    #[test]
    fn test_fitted_coeffs_multiple_temperatures() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();
        nasa.set_T_interval(900.0, 1100.0);
        nasa.norm_threshold = 1.0;

        assert!(nasa.fitting_coeffs_for_T_interval().is_ok());

        nasa.create_closures_Cp_dH_dS();
        nasa.create_sym_Cp_dH_dS();

        let test_temps = vec![920.0, 1000.0, 1080.0];

        for temp in test_temps {
            nasa.calculate_Cp_dH_dS(temp);

            // Test consistency between calculation methods
            let cp_closure = (nasa.C_fun)(temp);
            let dh_closure = (nasa.dh_fun)(temp);
            let ds_closure = (nasa.ds_fun)(temp);

            let cp_sym = nasa.Cp_sym.lambdify1D()(temp);
            let dh_sym = nasa.dh_sym.lambdify1D()(temp);
            let ds_sym = nasa.ds_sym.lambdify1D()(temp);

            assert_relative_eq!(cp_closure, nasa.Cp, epsilon = 1e-6);
            assert_relative_eq!(dh_closure, nasa.dh, epsilon = 1e-6);
            assert_relative_eq!(ds_closure, nasa.ds, epsilon = 1e-6);

            assert_relative_eq!(cp_sym, nasa.Cp, epsilon = 1e-6);
            assert_relative_eq!(dh_sym, nasa.dh, epsilon = 1e-6);
            assert_relative_eq!(ds_sym, nasa.ds, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_integr_mean_basic() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.extract_coefficients(500.0).unwrap();
        nasa.set_T_interval(400.0, 600.0);
        nasa.create_sym_Cp_dH_dS();
        match nasa.integr_mean() {
            Ok(_) => (),
            Err(e) => panic!("Error: {:?}", e),
        }
        // Calculate mean values
        assert!(nasa.integr_mean().is_ok());

        let mean_cp = nasa.Cp;
        let mean_dh = nasa.dh;
        let mean_ds = nasa.ds;

        // Verify mean values are reasonable
        assert!(mean_cp > 0.0);
        assert!(mean_dh != 0.0);
        assert!(mean_ds > 0.0);

        // Compare with point calculations at interval midpoint
        nasa.calculate_Cp_dH_dS(500.0);
        let midpoint_cp = nasa.Cp;
        let midpoint_dh = nasa.dh;
        let midpoint_ds = nasa.ds;

        // Mean values should be close to midpoint values for small intervals
        assert_relative_eq!(mean_cp, midpoint_cp, epsilon = 0.1);
        assert_relative_eq!(mean_dh, midpoint_dh, epsilon = 20_000.0);
        assert_relative_eq!(mean_ds, midpoint_ds, epsilon = 1000.0);
    }

    #[test]
    fn test_integr_mean_error_handling() {
        let mut nasa = NASAdata::new();

        // Test without T_interval set
        assert!(matches!(
            nasa.integr_mean(),
            Err(NASAError::InvalidTemperatureRange)
        ));

        // Test with invalid temperature range
        nasa.set_T_interval(600.0, 400.0); // T_min > T_max
        assert!(matches!(
            nasa.integr_mean(),
            Err(NASAError::InvalidTemperatureRange)
        ));

        // Test with equal temperatures
        nasa.set_T_interval(500.0, 500.0);
        assert!(matches!(
            nasa.integr_mean(),
            Err(NASAError::InvalidTemperatureRange)
        ));
    }

    #[test]
    fn test_fitting_non_adjacent() {
        let mut nasa = NASAdata::new();

        let coeffs = vec![
            Coeffs {
                T: (300.0, 800.0),
                coeff: (3.5, 0.001, -0.0001, 0.00001, -0.000001, -14000.0, 4.0),
            },
            Coeffs {
                T: (1200.0, 2000.0),
                coeff: (4.0, 0.0005, -0.00005, 0.000005, -0.0000005, -15000.0, 5.0),
            },
        ];

        nasa.norm_threshold = 1.0;
        let result = nasa.fitting_non_adjacent(coeffs, 500.0, 1500.0);
        assert!(result.is_ok());

        //    let report = result.unwrap();

        assert_ne!(nasa.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_fitting_non_adjacent_with_real_data() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();
        nasa.norm_threshold = 1.0;

        let coeffs_vec: Vec<_> = nasa.coeffs_map.values().cloned().collect();
        let result = nasa.fitting_non_adjacent(coeffs_vec, 400.0, 2500.0);

        assert!(result.is_ok());
        assert_ne!(nasa.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_fitting_non_adjacent_multiple_ranges() {
        let mut nasa = NASAdata::new();
        nasa.norm_threshold = 1.0;

        let coeffs = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: (3.0, 0.001, -0.0001, 0.00001, -0.000001, -14000.0, 4.0),
            },
            Coeffs {
                T: (800.0, 1200.0),
                coeff: (3.5, 0.0008, -0.00008, 0.000008, -0.0000008, -14500.0, 4.5),
            },
            Coeffs {
                T: (1500.0, 2000.0),
                coeff: (4.0, 0.0005, -0.00005, 0.000005, -0.0000005, -15000.0, 5.0),
            },
        ];

        let result = nasa.fitting_non_adjacent(coeffs, 400.0, 1800.0);
        assert!(result.is_ok());

        let report = result.unwrap();

        assert_ne!(nasa.coeffs, (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_fitting_non_adjacent_real_data_for_several_subs2() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let subs = vec!["H2O", "CO2", "NO", "N2"];
        for sub in subs {
            let sub_data = sublib.get(sub).unwrap();
            let mut nasa = NASAdata::new();
            nasa.from_serde(sub_data.clone()).unwrap();
            nasa.parse_coefficients().unwrap();
            println!("Sub: {}, number of T ranges {}", sub, nasa.coeffs_map.len());
        }
    }

    #[test]
    fn test_create_closures_with_T_range() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();

        // Test with temperature in first range
        assert!(nasa.create_closures_Cp_dH_dS_with_T_range(500.0).is_ok());
        let cp_val = (nasa.C_fun)(500.0);
        assert!(cp_val > 0.0);

        // Test with temperature in second range
        assert!(nasa.create_closures_Cp_dH_dS_with_T_range(1500.0).is_ok());
        let cp_val2 = (nasa.C_fun)(1500.0);
        assert!(cp_val2 > 0.0);
        assert_ne!(cp_val, cp_val2); // Should be different due to different coefficients
    }

    #[test]
    fn test_calculate_Cp_dH_dS_with_T_range() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();

        // Test with temperature in first range
        assert!(nasa.calculate_Cp_dH_dS_with_T_range(500.0).is_ok());
        let cp1 = nasa.Cp;
        let dh1 = nasa.dh;
        let ds1 = nasa.ds;
        assert!(cp1 > 0.0);

        // Test with temperature in second range
        assert!(nasa.calculate_Cp_dH_dS_with_T_range(1500.0).is_ok());
        let cp2 = nasa.Cp;
        let dh2 = nasa.dh;
        let ds2 = nasa.ds;
        assert!(cp2 > 0.0);

        // Values should be different due to different temperature ranges
        assert_ne!(cp1, cp2);
        assert_ne!(dh1, dh2);
        assert_ne!(ds1, ds2);
    }

    #[test]
    fn test_create_sym_Cp_dH_dS_with_T_range() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();

        // Test with temperature in first range
        assert!(nasa.create_sym_Cp_dH_dS_with_T_range(500.0).is_ok());
        let cp_sym1 = nasa.Cp_sym.clone();
        let cp_func1 = cp_sym1.lambdify1D();
        let cp_val1 = cp_func1(500.0);
        assert!(cp_val1 > 0.0);

        // Test with temperature in second range
        assert!(nasa.create_sym_Cp_dH_dS_with_T_range(1500.0).is_ok());
        let cp_sym2 = nasa.Cp_sym.clone();
        let cp_func2 = cp_sym2.lambdify1D();
        let cp_val2 = cp_func2(1500.0);
        assert!(cp_val2 > 0.0);

        // Symbolic expressions should be different due to different coefficients
        assert_ne!(cp_sym1.to_string(), cp_sym2.to_string());
    }

    #[test]
    fn test_T_range_methods_consistency() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("NASA_gas").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let mut nasa = NASAdata::new();

        nasa.from_serde(CO_data.clone()).unwrap();
        nasa.parse_coefficients().unwrap();

        let test_temp = 800.0;

        // Calculate using T_range method
        assert!(nasa.calculate_Cp_dH_dS_with_T_range(test_temp).is_ok());
        let cp_calc = nasa.Cp;

        // Create closures using T_range method
        assert!(
            nasa.create_closures_Cp_dH_dS_with_T_range(test_temp)
                .is_ok()
        );
        let cp_closure = (nasa.C_fun)(test_temp);

        // Create symbolic using T_range method
        assert!(nasa.create_sym_Cp_dH_dS_with_T_range(test_temp).is_ok());
        let cp_sym = nasa.Cp_sym.lambdify1D()(test_temp);

        // All methods should give consistent results
        assert_relative_eq!(cp_calc, cp_closure, epsilon = 1e-6);
        assert_relative_eq!(cp_calc, cp_sym, epsilon = 1e-6);
    }

    #[test]
    fn test_T_range_methods_error_handling() {
        let mut nasa = NASAdata::new();
        nasa.input.Cp = vec![300.0, 1000.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        nasa.parse_coefficients().unwrap();

        // Test with temperature outside all ranges
        assert!(nasa.create_closures_Cp_dH_dS_with_T_range(2000.0).is_err());
        assert!(nasa.calculate_Cp_dH_dS_with_T_range(2000.0).is_err());
        assert!(nasa.create_sym_Cp_dH_dS_with_T_range(2000.0).is_err());
    }
}
