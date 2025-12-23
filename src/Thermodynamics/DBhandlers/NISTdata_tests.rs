#[cfg(test)]
mod tests {

    use crate::Thermodynamics::DBhandlers::NIST_parser::{Phase, SearchType};
    use crate::Thermodynamics::DBhandlers::NISTdata::NISTdata;
    use crate::Thermodynamics::DBhandlers::thermo_api::{ThermoCalculator, create_thermal_by_name};
    use approx::assert_relative_eq;
    use std::{thread, time};
    #[allow(non_upper_case_globals)]
    const smtime: time::Duration = time::Duration::from_secs(5);
    use thread::sleep;
    #[test]
    fn test_thermo_calculator_nist() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        // Test newinstance
        //   assert!(nist.newinstance().is_ok());

        // Test from_serde

        // Test set_unit
        //   assert!(nist.set_unit(EnergyUnit::J).is_ok()) ;
        //   assert!(nist.set_unit(Some(EnergyUnit::Cal)).is_ok());

        // Test extract_model_coefficients
        // assert!(nist.extract_model_coefficients(400.0).is_ok());

        // Test calculate_Cp_dH_dS
        let r = nist.extract_coefficients(400.0);
        assert!(r.is_ok());
        let r = nist.calculate_Cp_dH_dS(400.0);
        println!("Cp: {}, dh: {}, ds: {}", nist.Cp, nist.dh, nist.ds);
        assert!(r.is_ok());
        assert!(nist.Cp > 0.0);
        assert!(nist.dh != 0.0);
        assert!(nist.ds != 0.0);

        // Test create_closures_Cp_dH_dS
        let r = nist.create_closures_Cp_dH_dS();
        assert!(r.is_ok());
        let t = 400.0;
        assert_relative_eq!((nist.C_fun)(t), nist.Cp, epsilon = 1e-6);
        assert_relative_eq!((nist.dh_fun)(t), nist.dh, epsilon = 1e-6);
        assert_relative_eq!((nist.ds_fun)(t), nist.ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        let r = nist.create_sym_Cp_dH_dS();
        assert!(r.is_ok());
        let Cp_sym = &nist.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, nist.Cp, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_cp_dh_ds(400.0, 3).unwrap();
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());
        nist.pretty_print_data().unwrap();
        // Test pretty_print_data
        assert!(nist.pretty_print_data().is_ok());
        sleep(smtime)
    }
    #[test]
    fn test_thermo_calculator_nist_simple_substance() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("O2".to_owned(), SearchType::All, Phase::Gas);

        // Test newinstance
        //   assert!(nist.newinstance().is_ok());

        // Test from_serde

        // Test set_unit
        //   assert!(nist.set_unit(EnergyUnit::J).is_ok()) ;
        //   assert!(nist.set_unit(Some(EnergyUnit::Cal)).is_ok());

        // Test extract_model_coefficients
        // assert!(nist.extract_model_coefficients(400.0).is_ok());

        // Test calculate_Cp_dH_dS
        let r = nist.extract_coefficients(400.0);
        assert!(r.is_ok());
        let r = nist.calculate_Cp_dH_dS(400.0);
        println!("Cp: {}, dh: {}, ds: {}", nist.Cp, nist.dh, nist.ds);
        assert!(r.is_ok());
        assert!(nist.Cp > 0.0);
        assert!(nist.dh != 0.0);
        assert!(nist.ds != 0.0);

        // Test create_closures_Cp_dH_dS
        let r = nist.create_closures_Cp_dH_dS();
        assert!(r.is_ok());
        let t = 400.0;
        assert_relative_eq!((nist.C_fun)(t), nist.Cp, epsilon = 1e-6);
        assert_relative_eq!((nist.dh_fun)(t), nist.dh, epsilon = 1e-6);
        assert_relative_eq!((nist.ds_fun)(t), nist.ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        let r = nist.create_sym_Cp_dH_dS();
        assert!(r.is_ok());
        let Cp_sym = &nist.Cp_sym;
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, nist.Cp, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_cp_dh_ds(400.0, 3).unwrap();
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());
        nist.pretty_print_data().unwrap();
        // Test pretty_print_data
        assert!(nist.pretty_print_data().is_ok());
        sleep(smtime)
    }

    #[test]
    fn test_thermo_calculator_nist_error_handling() {
        let mut nist = NISTdata::new();

        // Test invalid serde data - use data that would cause type mismatch
        let invalid_data = serde_json::json!({
            "cp": "not_a_vector",  // This should be Vec<Vec<f64>> but we're giving it a string
            "T": 123  // This should be Vec<Vec<f64>> but we're giving it a number
        });
        let result = nist.from_serde(invalid_data);

        assert!(result.is_err());
        sleep(smtime)
    }

    #[test]
    fn test_nist_clone() {
        let mut nist = NISTdata::new();

        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        let _ = nist.extract_coefficients(400.0);
        let _ = nist.calculate_Cp_dH_dS(400.0);
        let _ = nist.create_closures_Cp_dH_dS();
        let _ = nist.create_sym_Cp_dH_dS();

        // Clone the instance
        let mut nist_clone = nist.clone();
        let _ = nist_clone.extract_coefficients(400.0);
        // Test that the clone has the same values
        assert_relative_eq!(nist_clone.Cp, nist.Cp, epsilon = 1e-6);
        assert_relative_eq!(nist_clone.dh, nist.dh, epsilon = 1e-6);
        assert_relative_eq!(nist_clone.ds, nist.ds, epsilon = 1e-6);

        // Test that the cloned functions work
        let t = 400.0;
        assert_relative_eq!((nist_clone.C_fun)(t), (nist.C_fun)(t), epsilon = 1e-6);
        assert_relative_eq!((nist_clone.dh_fun)(t), (nist.dh_fun)(t), epsilon = 1e-6);
        assert_relative_eq!((nist_clone.ds_fun)(t), (nist.ds_fun)(t), epsilon = 1e-6);

        // Test that the symbolic expressions are the same
        assert_eq!(nist_clone.Cp_sym.to_string(), nist.Cp_sym.to_string());
        assert_eq!(nist_clone.dh_sym.to_string(), nist.dh_sym.to_string());
        assert_eq!(nist_clone.ds_sym.to_string(), nist.ds_sym.to_string());
        sleep(smtime)
    }
    #[test]
    fn test_temperature_range_calculations() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CH4".to_owned(), SearchType::All, Phase::Gas);

        // Test calculations at different temperatures
        let temperatures = [298.15, 400.0, 600.0, 800.0, 1000.0];
        for &T in &temperatures {
            if nist.extract_coefficients(T).is_ok() {
                let _ = nist.calculate_cp_dh_ds(T);
                assert!(nist.Cp > 0.0, "Cp should be positive at {} K", T);
                assert!(nist.dh != 0.0, "dh should not be zero at {} K", T);
                assert!(nist.ds != 0.0, "ds should not be zero at {} K", T);
            }
        }
        sleep(smtime)
    }

    #[test]
    fn test_unit_conversion() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);
        let _ = nist.extract_coefficients(400.0);

        // Test Joules (default)
        let _ = nist.set_unit("J");
        let _ = nist.calculate_cp_dh_ds(400.0);
        let cp_j = nist.Cp;

        // Test Calories
        let _ = nist.set_unit("cal");
        let _ = nist.calculate_cp_dh_ds(400.0);
        let cp_cal = nist.Cp;

        // Verify conversion (1 cal = 4.184 J)
        assert_relative_eq!(cp_j / cp_cal, 4.184, epsilon = 0.01);
        sleep(smtime)
    }

    #[test]
    fn test_symbolic_expressions() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);
        let _ = nist.extract_coefficients(500.0);
        let _ = nist.create_sym_cp_dh_ds();

        // Test symbolic expressions are not empty
        assert!(!nist.Cp_sym.is_zero());
        assert!(!nist.dh_sym.is_zero());
        assert!(!nist.ds_sym.is_zero());

        // Test symbolic evaluation matches numerical
        let _ = nist.calculate_cp_dh_ds(500.0);
        let cp_num = nist.Cp;
        let cp_sym_eval = nist.Cp_sym.lambdify1D()(500.0);
        assert_relative_eq!(cp_num, cp_sym_eval, epsilon = 1e-6);
        sleep(smtime)
    }

    #[test]
    fn test_closure_functions() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("NO2".to_owned(), SearchType::All, Phase::Gas);
        let _ = nist.extract_coefficients(350.0);
        let _ = nist.create_closure_cp_dh_ds();

        // Test closures at multiple temperatures
        let test_temps = [300.0, 350.0, 400.0, 450.0];
        for &T in &test_temps {
            let cp_closure = (nist.C_fun)(T);
            let dh_closure = (nist.dh_fun)(T);
            let ds_closure = (nist.ds_fun)(T);

            assert!(cp_closure > 0.0, "Cp closure should be positive at {} K", T);
            assert!(
                dh_closure != 0.0,
                "dh closure should not be zero at {} K",
                T
            );
            assert!(
                ds_closure != 0.0,
                "ds closure should not be zero at {} K",
                T
            );
        }
        sleep(smtime)
    }

    #[test]
    fn test_taylor_series_expansion() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("NO2".to_owned(), SearchType::All, Phase::Gas);
        let _ = nist.extract_coefficients(400.0);
        let _ = nist.create_sym_cp_dh_ds();

        // Test Taylor series at different orders
        for order in 1..=5 {
            let result = nist.clone().Taylor_series_Cp_dH_dS(400.0, order);
            assert!(
                result.is_ok(),
                "Taylor series should work for order {}",
                order
            );

            let (cp_taylor, dh_taylor, ds_taylor) = result.unwrap();
            assert!(!cp_taylor.is_zero(), "Cp Taylor series should not be zero");
            assert!(!dh_taylor.is_zero(), "dh Taylor series should not be zero");
            assert!(!ds_taylor.is_zero(), "ds Taylor series should not be zero");
        }
        sleep(smtime)
    }

    #[test]
    fn test_coefficient_extraction() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        // Test valid temperature
        assert!(nist.extract_coefficients(400.0).is_ok());
        assert!(nist.coeffs.is_some());

        let _coeffs = nist.coeffs.unwrap();
        // assert_eq!(coeffs.len(), 8); // Should have 8 coefficients (a,b,c,d,e,f,g,h)

        // Test invalid temperature (outside range)
        let result = nist.extract_coefficients(10000.0);
        assert!(result.is_err());
        sleep(smtime)
    }

    #[test]
    fn test_integration_mean() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);
        let _ = nist.extract_coefficients(500.0);
        let _ = nist.create_sym_cp_dh_ds();

        // Set temperature interval
        nist.set_T_interval(400.0, 600.0);

        // Test integration
        match nist.integr_mean() {
            Ok(_) => {}
            Err(e) => panic!("Integration failed: {}", e),
        }

        // Check that mean values are calculated
        assert!(nist.Cp > 0.0);
        assert!(nist.dh != 0.0);
        assert!(nist.ds != 0.0);
        sleep(smtime)
    }

    #[test]
    fn ThermoCalculator_nist() {
        let mut nist = create_thermal_by_name("NIST");
        let _ = nist.newinstance();
        let _ = nist.renew_base("CO".to_owned(), SearchType::All, Phase::Gas);
        let T = 400.0;
        let _ = nist.extract_model_coefficients(T);
        let _ = nist.calculate_Cp_dH_dS(400.0);
        let Cp = nist.get_Cp().unwrap();
        let dh = nist.get_dh().unwrap();
        let ds = nist.get_ds().unwrap();
        assert!(Cp > 0.0);
        assert!(dh != 0.0);
        assert!(ds != 0.0);

        // Test create_closures_Cp_dH_dS
        let _ = nist.create_closures_Cp_dH_dS();
        let t = 400.0;
        assert_relative_eq!((nist.get_C_fun().unwrap())(t), Cp, epsilon = 1e-6);
        assert_relative_eq!((nist.get_dh_fun().unwrap())(t), dh, epsilon = 1e-6);
        assert_relative_eq!((nist.get_ds_fun().unwrap())(t), ds, epsilon = 1e-6);

        // Test create_sym_Cp_dH_dS
        let _ = nist.create_sym_Cp_dH_dS();
        let Cp_sym = &nist.get_Cp_sym().unwrap();
        let Cp_T = Cp_sym.lambdify1D();
        let Cp_value = Cp_T(400.0);
        assert_relative_eq!(Cp_value, Cp, epsilon = 1e-6);

        // Test Taylor_series_cp_dh_ds
        let (Cp_taylor, dh_taylor, ds_taylor) = nist.Taylor_series_cp_dh_ds(400.0, 3).unwrap();
        assert!(!Cp_taylor.is_zero());
        assert!(!dh_taylor.is_zero());
        assert!(!ds_taylor.is_zero());

        // Test pretty_print_data
        assert!(nist.pretty_print_data().is_ok());
        sleep(smtime)
    }

    ////////////////////////////TESTS FOR DATA FITTING//////////////////////////////////////////////

    #[test]
    fn test_fitting_coeffs_for_T_interval_same_interval() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();
        nist.set_T_interval(400.0, 600.0);

        assert!(nist.fitting_coeffs_for_T_interval().is_ok());
        assert!(nist.coeffs.is_some());
        sleep(smtime)
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_adjacent() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let first_range = nist.coeffs_map.get(&0).unwrap().T;
            let second_range = nist.coeffs_map.get(&1).unwrap().T;

            if first_range.1 == second_range.0 {
                let T_center = first_range.1;
                let T_min = T_center - 100.0;
                let T_max = T_center + 100.0;

                let _ = nist.extract_coefficients(T_min);
                let _ = nist.create_closure_cp_dh_ds();
                let mut T = T_min;
                let step = 10.0;
                let mut Cp: Vec<f64> = Vec::new();
                let mut dh: Vec<f64> = Vec::new();
                let mut ds: Vec<f64> = Vec::new();

                while T <= T_center {
                    let Cpi = (nist.C_fun)(T);
                    Cp.push(Cpi);
                    let dhi = (nist.dh_fun)(T);
                    dh.push(dhi);
                    let dsi = (nist.ds_fun)(T);
                    ds.push(dsi);
                    T += step;
                }
                let mut T = T_center + 10.0;
                let _ = nist.extract_coefficients(T_max);
                let _ = nist.create_closure_cp_dh_ds();
                while T <= T_max {
                    let Cpi = (nist.C_fun)(T);
                    Cp.push(Cpi);
                    let dhi = (nist.dh_fun)(T);
                    dh.push(dhi);
                    let dsi = (nist.ds_fun)(T);
                    ds.push(dsi);
                    T += step;
                }
                nist.set_T_interval(T_min, T_max);

                let result = nist.fitting_coeffs_for_T_interval();
                assert!(result.is_ok());
                assert!(nist.coeffs.is_some());
                let _ = nist.create_closure_cp_dh_ds();
                let mut T = T_min;

                let mut i = 0;
                while T <= T_max {
                    let Cp_fitted = (nist.C_fun)(T);
                    let delta_Cp = (Cp_fitted - Cp[i]) / Cp_fitted;

                    assert!(
                        delta_Cp.abs() < 0.1,
                        "Cp_fitted {}, Cp, {}",
                        Cp_fitted,
                        Cp[i]
                    );
                    let dh_fitted = (nist.dh_fun)(T);
                    let delta_dh = (dh_fitted - dh[i]) / dh_fitted;
                    assert!(delta_dh.abs() < 0.1);

                    let ds_fitted = (nist.ds_fun)(T);
                    let delta_ds = (ds_fitted - ds[i]) / ds_fitted;
                    assert!(delta_ds.abs() < 0.1);
                    T += step;
                    i += 1;
                }
            }
        }
        sleep(smtime)
    }

    #[test]
    fn test_fitting_coeffs_for_T_interval_non_adjacent() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 3 {
            let first_range = nist.coeffs_map.get(&0).unwrap().T;
            let third_range = nist.coeffs_map.get(&2).unwrap().T;

            let T_min = first_range.0 + 50.0;
            let T_max = third_range.1 - 50.0;
            nist.set_T_interval(T_min, T_max);

            let result = nist.fitting_coeffs_for_T_interval();
            assert!(result.is_ok());
            assert!(nist.coeffs.is_some());
        }
        sleep(smtime)
    }

    #[test]
    fn test_interval_for_this_T() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("NO2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if !nist.coeffs_map.is_empty() {
            let first_range = nist.coeffs_map.get(&0).unwrap().T;
            let test_temp = (first_range.0 + first_range.1) / 2.0;

            let result = nist.interval_for_this_T(test_temp);
            assert!(result.is_ok());

            let (interval_idx, coeffs) = result.unwrap();
            assert_eq!(interval_idx, 0);
            assert!(coeffs.T.0 <= test_temp && test_temp <= coeffs.T.1);
        }
        sleep(smtime)
    }

    #[test]
    fn test_fitted_coeffs_consistency() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("O2".to_owned(), SearchType::All, Phase::Gas);
        nist.parse_coefficients().unwrap();

        if nist.coeffs_map.len() >= 2 {
            let first_range = nist.coeffs_map.get(&0).unwrap().T;
            let second_range = nist.coeffs_map.get(&1).unwrap().T;

            if first_range.1 == second_range.0 {
                let T_center = first_range.1;
                let T_min = T_center - 50.0;
                let T_max = T_center + 50.0;
                nist.set_T_interval(T_min, T_max);

                assert!(nist.fitting_coeffs_for_T_interval().is_ok());

                let _ = nist.create_closure_cp_dh_ds();
                let _ = nist.create_sym_cp_dh_ds();

                let test_temp = T_center;
                let _ = nist.calculate_cp_dh_ds(test_temp);

                let cp_calc = nist.Cp;
                let dh_calc = nist.dh;
                let ds_calc = nist.ds;

                let cp_closure = (nist.C_fun)(test_temp);
                let dh_closure = (nist.dh_fun)(test_temp);
                let ds_closure = (nist.ds_fun)(test_temp);

                assert_relative_eq!(cp_closure, cp_calc, epsilon = 1e-6);
                assert_relative_eq!(dh_closure, dh_calc, epsilon = 1e-6);
                assert_relative_eq!(ds_closure, ds_calc, epsilon = 1e-6);
            }
        }
        sleep(smtime)
    }

    #[test]
    fn test_fitting_error_handling() {
        let mut nist = NISTdata::new();

        let result = nist.fitting_coeffs_for_T_interval();
        assert!(result.is_err());

        nist.set_T_interval(1000.0, 500.0);
        let result = nist.fitting_coeffs_for_T_interval();
        assert!(result.is_err());
        sleep(smtime)
    }

    #[test]
    fn test_fitting_non_adjacent_identical_coeffs() {
        use crate::Thermodynamics::DBhandlers::NISTdata::Coeffs;

        let mut nist = NISTdata::new();
        let identical_coeffs = (30.0, -3.0, 0.5, -0.05, 0.005, -800.0, 180.0, 0.0);

        let coeffs = vec![
            Coeffs {
                T: (300.0, 600.0),
                coeff: identical_coeffs,
            },
            Coeffs {
                T: (800.0, 1100.0),
                coeff: identical_coeffs,
            },
            Coeffs {
                T: (1300.0, 1600.0),
                coeff: identical_coeffs,
            },
        ];

        let result = nist.fitting_non_adjacent(coeffs, 400.0, 1400.0);
        assert!(result.is_ok());

        let fitted_coeffs = nist.coeffs.unwrap();
        assert_relative_eq!(fitted_coeffs.0, identical_coeffs.0, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.1, identical_coeffs.1, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.2, identical_coeffs.2, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.3, identical_coeffs.3, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.4, identical_coeffs.4, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.5, identical_coeffs.5, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.6, identical_coeffs.6, epsilon = 1e-3);
        assert_relative_eq!(fitted_coeffs.7, identical_coeffs.7, epsilon = 1e-3);
        sleep(smtime)
    }

    #[test]
    fn test_create_closures_with_T_range() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO".to_owned(), SearchType::All, Phase::Gas);

        // Test with temperature in first range
        assert!(nist.create_closures_Cp_dH_dS_with_T_range(400.0).is_ok());
        let cp_val = (nist.C_fun)(400.0);
        assert!(cp_val > 0.0);

        // Test with temperature in different range if available
        if nist.coeffs_map.len() > 1 {
            let second_range = nist.coeffs_map.get(&1).unwrap().T;
            let test_temp = (second_range.0 + second_range.1) / 2.0;
            assert!(
                nist.create_closures_Cp_dH_dS_with_T_range(test_temp)
                    .is_ok()
            );
            let cp_val2 = (nist.C_fun)(test_temp);
            assert!(cp_val2 > 0.0);
        }
        sleep(smtime)
    }

    #[test]
    fn test_calculate_Cp_dH_dS_with_T_range() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("H2O".to_owned(), SearchType::All, Phase::Gas);

        // Test with temperature in first range
        assert!(nist.calculate_Cp_dH_dS_with_T_range(500.0).is_ok());
        let cp1 = nist.Cp;
        let dh1 = nist.dh;
        let ds1 = nist.ds;
        assert!(cp1 > 0.0);

        // Test with temperature in different range if available
        if nist.coeffs_map.len() > 1 {
            let second_range = nist.coeffs_map.get(&1).unwrap().T;
            let test_temp = (second_range.0 + second_range.1) / 2.0;
            assert!(nist.calculate_Cp_dH_dS_with_T_range(test_temp).is_ok());
            let cp2 = nist.Cp;
            let dh2 = nist.dh;
            let ds2 = nist.ds;
            assert!(cp2 > 0.0);

            // Values should be different due to different temperature ranges
            assert_ne!(cp1, cp2);
            assert_ne!(dh1, dh2);
            assert_ne!(ds1, ds2);
        }
        sleep(smtime)
    }

    #[test]
    fn test_create_sym_Cp_dH_dS_with_T_range() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CO2".to_owned(), SearchType::All, Phase::Gas);

        // Test with temperature in first range
        assert!(nist.create_sym_Cp_dH_dS_with_T_range(400.0).is_ok());
        let cp_sym1 = nist.Cp_sym.clone();
        let cp_func1 = cp_sym1.lambdify1D();
        let cp_val1 = cp_func1(400.0);
        assert!(cp_val1 > 0.0);

        // Test with temperature in different range if available
        if nist.coeffs_map.len() > 1 {
            let second_range = nist.coeffs_map.get(&1).unwrap().T;
            let test_temp = (second_range.0 + second_range.1) / 2.0;
            assert!(nist.create_sym_Cp_dH_dS_with_T_range(test_temp).is_ok());
            let cp_sym2 = nist.Cp_sym.clone();
            let cp_func2 = cp_sym2.lambdify1D();
            let cp_val2 = cp_func2(test_temp);
            assert!(cp_val2 > 0.0);

            // Symbolic expressions should be different due to different coefficients
            assert_ne!(cp_sym1.to_string(), cp_sym2.to_string());
        }
        sleep(smtime)
    }

    #[test]
    fn test_T_range_methods_consistency() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("NO2".to_owned(), SearchType::All, Phase::Gas);

        let test_temp = 400.0;

        // Calculate using T_range method
        assert!(nist.calculate_Cp_dH_dS_with_T_range(test_temp).is_ok());
        let cp_calc = nist.Cp;

        // Create closures using T_range method
        assert!(
            nist.create_closures_Cp_dH_dS_with_T_range(test_temp)
                .is_ok()
        );
        let cp_closure = (nist.C_fun)(test_temp);

        // Create symbolic using T_range method
        assert!(nist.create_sym_Cp_dH_dS_with_T_range(test_temp).is_ok());
        let cp_sym = nist.Cp_sym.lambdify1D()(test_temp);

        // All methods should give consistent results
        assert_relative_eq!(cp_calc, cp_closure, epsilon = 1e-6);
        assert_relative_eq!(cp_calc, cp_sym, epsilon = 1e-6);
        sleep(smtime)
    }

    #[test]
    fn test_T_range_methods_error_handling() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("O2".to_owned(), SearchType::All, Phase::Gas);

        // Test with temperature outside all ranges
        let result1 = nist.create_closures_Cp_dH_dS_with_T_range(10000.0);
        let result2 = nist.calculate_Cp_dH_dS_with_T_range(10000.0);
        let result3 = nist.create_sym_Cp_dH_dS_with_T_range(10000.0);

        assert!(result1.is_err());
        assert!(result2.is_err());
        assert!(result3.is_err());
        sleep(smtime)
    }

    #[test]
    fn test_T_range_methods_multiple_calls() {
        let mut nist = NISTdata::new();
        let _ = nist.get_data_from_NIST("CH4".to_owned(), SearchType::All, Phase::Gas);

        let temps = [350.0, 450.0, 550.0];

        for &temp in &temps {
            // Test that multiple calls work correctly
            assert!(nist.calculate_Cp_dH_dS_with_T_range(temp).is_ok());
            assert!(nist.create_closures_Cp_dH_dS_with_T_range(temp).is_ok());
            assert!(nist.create_sym_Cp_dH_dS_with_T_range(temp).is_ok());

            // Verify consistency between methods
            let cp_calc = nist.Cp;
            let cp_closure = (nist.C_fun)(temp);
            let cp_sym = nist.Cp_sym.lambdify1D()(temp);

            assert_relative_eq!(cp_calc, cp_closure, epsilon = 1e-6);
            assert_relative_eq!(cp_calc, cp_sym, epsilon = 1e-6);
        }
        sleep(smtime)
    }
}
