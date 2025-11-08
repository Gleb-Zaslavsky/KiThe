#[cfg(test)]
mod tests {
    use crate::Thermodynamics::DBhandlers::NIST_parser::{
        NistError, NistInput, NistParser, Phase, SearchType,
    };
    // use std::collections::HashMap;
    /*
       // Mock HTTP client for testing
       #[derive(Default)]
       struct MockHttpClient {
           responses: HashMap<String, Result<String, reqwest::Error>>,
       }

       impl MockHttpClient {
           fn new() -> Self {
               Self {
                   responses: HashMap::new(),
               }
           }

           fn mock_response(&mut self, url: &str, response: Result<String, reqwest::Error>) {
               self.responses.insert(url.to_string(), response);
           }
       }

       impl HttpClient for MockHttpClient {
           fn get_text(&self, url: &str) -> Result<String, reqwest::Error> {
               let resp = self.responses.get(url).unwrap().clone()?;
               match resp {
                   Ok(text) => Ok(text.clone()),
                   Err(err) => Err(err),
               }
           }
       }
    */
    #[test]
    fn test_url_construction() {
        let parser = NistParser::new();

        // Test CAS number
        let cas_url = parser.construct_url("7732-18-5").unwrap();
        assert_eq!(
            cas_url.as_str(),
            "https://webbook.nist.gov/cgi/cbook.cgi?ID=7732-18-5&Units=SI"
        );

        // Test chemical formula
        let formula_url = parser.construct_url("H2O").unwrap();
        assert_eq!(
            formula_url.as_str(),
            "https://webbook.nist.gov/cgi/cbook.cgi?Formula=H2O&NoIon=on&Units=SI"
        );

        // Test chemical name
        let name_url = parser.construct_url("Water").unwrap();
        assert_eq!(
            name_url.as_str(),
            "https://webbook.nist.gov/cgi/cbook.cgi?Name=Water&Units=SI"
        );

        // Test space handling
        let space_url = parser.construct_url("Carbon dioxide").unwrap();
        assert_eq!(
            space_url.as_str(),
            "https://webbook.nist.gov/cgi/cbook.cgi?Name=Carbondioxide&Units=SI"
        );
    }

    #[test]
    fn test_phase_conversion() {
        assert_eq!(Phase::Gas.as_str(), "gas");
        assert_eq!(Phase::Solid.as_str(), "solid");
    }

    #[test]
    fn test_real_substance_fetch() {
        let parser = NistParser::new();

        // Test methane data fetch
        let result = parser.get_data("CH4", SearchType::MolarMass, Phase::Gas);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.molar_mass.is_some());
            assert!(data.cp.is_none()); // We didn't request Cp
            assert!(data.dh.is_none()); // We didn't request dH
        }

        // Test water data fetch
        let result = parser.get_data("H2O", SearchType::Cp, Phase::Gas);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.cp.is_some());
            assert!(data.molar_mass.is_none()); // We didn't request molar mass
        }
    }

    #[test]
    fn test_nonexistent_substance() {
        let parser = NistParser::new();
        let result = parser.get_data("NonexistentSubstance123", SearchType::MolarMass, Phase::Gas);
        assert!(matches!(result, Err(NistError::SubstanceNotFound)));
    }

    #[test]
    fn test_thermodynamic_data() {
        let parser = NistParser::new();

        // Test enthalpy data for water
        let result = parser.get_data("H2O", SearchType::DeltaH, Phase::Gas);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.dh.is_some() || data.ds.is_some());
        }
    }

    #[test]
    fn test_data_structure_initialization() {
        let data = NistInput {
            cp: Some(vec![vec![1.0, 2.0, 3.0]]),
            T: Some(vec![vec![298.0, 1000.0]]),
            dh: Some(42.0),
            ds: Some(100.0),
            coeffs: None,
            molar_mass: Some(18.015),
            unit: None,
            unit_multiplier: 1.0,
        };

        assert_eq!(data.cp.as_ref().unwrap()[0][0], 1.0);
        assert_eq!(data.dh.unwrap(), 42.0);
        assert_eq!(data.ds.unwrap(), 100.0);
        assert_eq!(data.molar_mass.unwrap(), 18.015);
    }

    #[test]
    fn test_ch4_gas() {
        let parser = NistParser::new();
        let substance = "CH4";
        let T = 300.0;
        let mut result = parser.get_data(substance, SearchType::All, Phase::Gas);
        assert!(result.is_ok());
        if let Ok(data) = result.as_mut() {
            let _ = data.extract_coefficients(T);
            let (Cp, dh, ds) = data
                .caclc_cp_dh_ds(T)
                .expect("Error calculating cp, dh, ds");
            assert!(Cp > 0.0);
            assert!(dh < 0.0); // Assuming formation enthalpy is negative
            assert!(ds > 0.0);
            println!(
                "CH4 Gas: Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                Cp, dh, ds
            );
        }
    }

    #[test]
    fn test_nacl_solid() {
        let parser = NistParser::new();
        let substance = "NaCl";
        let T = 300.0;
        let mut result = parser.get_data(substance, SearchType::All, Phase::Solid);
        assert!(result.is_ok());
        if let Ok(data) = result.as_mut() {
            let _ = data.extract_coefficients(T);
            let (Cp, dh, ds) = data
                .caclc_cp_dh_ds(T)
                .expect("Error calculating cp, dh, ds");
            assert!(Cp > 0.0);
            assert!(dh < 0.0); // Assuming formation enthalpy is negative
            assert!(ds > 0.0);
            println!(
                "NaCl Solid: Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                Cp, dh, ds
            );
        }
    }

    #[test]
    fn test_nacl_liquid() {
        let parser = NistParser::new();
        let substance = "NaCl";
        let T = 1200.15;
        let mut result = parser.get_data(substance, SearchType::All, Phase::Liquid);
        assert!(result.is_ok());

        if let Ok(data) = result.as_mut() {
            let _ = data.extract_coefficients(T);
            let (Cp, dh, ds) = data
                .caclc_cp_dh_ds(1200.15)
                .expect("Error calculating cp, dh, ds");
            assert!(Cp > 0.0);
            assert!(dh < 0.0); // Assuming formation enthalpy is negative
            assert!(ds > 0.0);
            println!(
                "NaCl Liquid: Cp J/mol*K: {}, dh kJ/mol: {}, ds J/mol*K: {}",
                Cp, dh, ds
            );
        }
    }

    #[test]
    fn test_closures() {
        let parser = NistParser::new();
        let substance = "H2O";
        let T = 300.0;
        let mut result = parser.get_data(substance, SearchType::All, Phase::Liquid);
        assert!(result.is_ok());
        if let Ok(data) = result.as_mut() {
            let _ = data.extract_coefficients(T);
            let (Cp, dh, ds) = data
                .caclc_cp_dh_ds(T)
                .expect("Error calculating cp, dh, ds");
            let (Cp_fun, dh_fun, ds_fun) = data
                .create_closure_cp_dh_ds()
                .expect("Error calculating cp, dh, ds");
            let Cp_eval = Cp_fun(T);
            let dh_eval = dh_fun(T);
            let ds_eval = ds_fun(T);
            println!(
                " Cp J/mol*K: {}, {}, dh kJ/mol: {}, {} ds J/mol*K: {} {}",
                Cp, Cp_eval, dh, dh_eval, ds, ds_eval
            );
            assert_eq!(Cp_eval, Cp);
            assert_eq!(dh_eval, dh);
            assert_eq!(ds_eval, ds);
        }
    }

    #[test]
    fn test_sym_functions() {
        let parser = NistParser::new();
        let substance = "H2O";
        let T = 300.0;
        let mut result = parser.get_data(substance, SearchType::All, Phase::Liquid);
        assert!(result.is_ok());
        if let Ok(data) = result.as_mut() {
            let _ = data.extract_coefficients(T);
            let (Cp, dh, ds) = data
                .caclc_cp_dh_ds(T)
                .expect("Error calculating cp, dh, ds");
            let (Cp_sym, dh_sym, ds_sym) = data
                .create_sym_cp_dh_ds()
                .expect("Error calculating cp, dh, ds");
            let Cp_eval = Cp_sym.lambdify1D()(T);
            let dh_eval = dh_sym.lambdify1D()(T);
            let ds_eval = ds_sym.lambdify1D()(T);
            println!(
                " Cp J/mol*K: {}, {}, dh kJ/mol: {}, {} ds J/mol*K: {} {}",
                Cp, Cp_eval, dh, dh_eval, ds, ds_eval
            );
            assert_eq!(Cp_eval, Cp);
            assert_eq!(dh_eval, dh);
            assert_eq!(ds_eval, ds);
        }
    }
    /*
    #[test]
    fn test_mock_substance_fetch() {
        let mut mock_client = MockHttpClient::new();

        // Mock response for methane
        let mock_response = r#"
            <html>
                <body>
                    <h1 id="Top">Methane</h1>
                    <li>Molecular weight: 16.043</li>
                    <table aria-label="Constant pressure heat capacity of gas">
                        <tr><td>298.15</td><td>35.69</td></tr>
                    </table>
                </body>
            </html>
        "#;

        mock_client.mock_response(
            "https://webbook.nist.gov/cgi/cbook.cgi?Formula=CH4&NoIon=on&Units=SI",
            Ok(mock_response.to_string())
        );

        let parser = NistParser::with_client(mock_client);
        let result = parser.get_data("CH4", SearchType::MolarMass, Phase::Gas);

        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.molar_mass.is_some());
            if let Some(mass) = data.molar_mass {
                assert!((mass - 16.043).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_mock_not_found() {
        let mut mock_client = MockHttpClient::new();

        // Mock response for non-existent substance
        let mock_response = r#"
            <html>
                <body>
                    <h1>Name Not Found</h1>
                </body>
            </html>
        "#;

        mock_client.mock_response(
            "https://webbook.nist.gov/cgi/cbook.cgi?Name=NonexistentSubstance&Units=SI",
            Ok(mock_response.to_string())
        );

        let parser = NistParser::with_client(mock_client);
        let result = parser.get_data("NonexistentSubstance", SearchType::MolarMass, Phase::Gas);

        assert!(matches!(result, Err(NistError::SubstanceNotFound)));
    }

    #[test]
    fn test_mock_network_error() {
        let mut mock_client = MockHttpClient::new();

        // Mock a network error
        mock_client.mock_response(
            "https://webbook.nist.gov/cgi/cbook.cgi?Formula=H2O&NoIon=on&Units=SI",
            Err(reqwest::Error::from(std::io::Error::new(
                std::io::ErrorKind::ConnectionRefused,
                "Connection refused"
            )))
        );

        let parser = NistParser::with_client(mock_client);
        let result = parser.get_data("H2O", SearchType::MolarMass, Phase::Gas);

        assert!(matches!(result, Err(NistError::NetworkError(_))));
    }
    */
}
