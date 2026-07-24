#[cfg(test)]
mod tests {
    use crate::Thermodynamics::DBhandlers::NIST_parser::{
        NistError, NistInput, NistParser, Phase, SearchType,
    };
    use reqwest::blocking::Client;
    use std::io::{Read, Write};
    use std::net::TcpListener;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };
    use std::thread;
    use std::time::Duration;
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

    fn build_timeout_client(timeout_ms: u64) -> Client {
        Client::builder()
            .connect_timeout(Duration::from_millis(timeout_ms))
            .timeout(Duration::from_millis(timeout_ms))
            .build()
            .expect("test client should be constructible")
    }

    fn http_response(status: &str, body: &str) -> String {
        format!(
            "HTTP/1.1 {}\r\nContent-Length: {}\r\nContent-Type: text/html; charset=utf-8\r\nConnection: close\r\n\r\n{}",
            status,
            body.len(),
            body
        )
    }

    fn spawn_response_server(response: String) -> (String, Arc<AtomicUsize>) {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind test server");
        let addr = listener.local_addr().expect("local addr");
        let hits = Arc::new(AtomicUsize::new(0));
        let hits_for_thread = Arc::clone(&hits);

        thread::spawn(move || {
            if let Ok((mut stream, _)) = listener.accept() {
                hits_for_thread.fetch_add(1, Ordering::SeqCst);
                let mut buf = [0_u8; 1024];
                let _ = stream.read(&mut buf);
                let _ = stream.write_all(response.as_bytes());
                let _ = stream.flush();
            }
        });

        (format!("http://{}", addr), hits)
    }

    fn spawn_hanging_server() -> (String, Arc<AtomicUsize>) {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind test server");
        let addr = listener.local_addr().expect("local addr");
        let hits = Arc::new(AtomicUsize::new(0));
        let hits_for_thread = Arc::clone(&hits);

        thread::spawn(move || {
            if let Ok((mut stream, _)) = listener.accept() {
                hits_for_thread.fetch_add(1, Ordering::SeqCst);
                let mut buf = [0_u8; 1024];
                let _ = stream.read(&mut buf);
                thread::sleep(Duration::from_secs(2));
            }
        });

        (format!("http://{}", addr), hits)
    }

    fn spawn_threshold_server(
        success_response: String,
        success_requests: usize,
    ) -> (String, Arc<AtomicUsize>) {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind test server");
        let addr = listener.local_addr().expect("local addr");
        let hits = Arc::new(AtomicUsize::new(0));
        let hits_for_thread = Arc::clone(&hits);
        let success_response = Arc::new(success_response);
        let success_response_for_thread = Arc::clone(&success_response);

        thread::spawn(move || {
            for _ in 0..32 {
                if let Ok((mut stream, _)) = listener.accept() {
                    let current = hits_for_thread.fetch_add(1, Ordering::SeqCst) + 1;
                    let mut buf = [0_u8; 1024];
                    let _ = stream.read(&mut buf);
                    let response = if current <= success_requests {
                        success_response_for_thread.as_str().to_string()
                    } else {
                        http_response("404 Not Found", "missing")
                    };
                    let _ = stream.write_all(response.as_bytes());
                    let _ = stream.flush();
                } else {
                    break;
                }
            }
        });

        (format!("http://{}", addr), hits)
    }

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
    fn test_fetch_page_rejects_http_error_status() {
        let (base_url, hits) = spawn_response_server(http_response("404 Not Found", "missing"));
        let client = build_timeout_client(250);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let url = url::Url::parse(&format!("{}/anything", base_url)).unwrap();

        let result = parser.fetch_page(&url);

        assert_eq!(hits.load(Ordering::SeqCst), 1);
        assert!(matches!(result, Err(NistError::HttpStatus(status)) if status.as_u16() == 404));
    }

    #[test]
    fn test_fetch_page_rejects_plain_text_body() {
        let (base_url, _) = spawn_response_server(http_response("200 OK", "plain text only"));
        let client = build_timeout_client(250);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let url = url::Url::parse(&format!("{}/anything", base_url)).unwrap();

        let result = parser.fetch_page(&url);

        assert!(matches!(result, Err(NistError::InvalidDataFormat)));
    }

    #[test]
    fn test_fetch_page_rejects_oversized_body() {
        let body = format!(
            "<html><body>{}</body></html>",
            "x".repeat(2 * 1024 * 1024 + 64)
        );
        let (base_url, _) = spawn_response_server(http_response("200 OK", &body));
        let client = build_timeout_client(250);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let url = url::Url::parse(&format!("{}/anything", base_url)).unwrap();

        let result = parser.fetch_page(&url);

        assert!(matches!(result, Err(NistError::ResponseTooLarge { .. })));
    }

    #[test]
    fn test_fetch_page_times_out_on_hanging_server() {
        let (base_url, hits) = spawn_hanging_server();
        let client = build_timeout_client(100);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let url = url::Url::parse(&format!("{}/anything", base_url)).unwrap();

        let result = parser.fetch_page(&url);

        assert_eq!(hits.load(Ordering::SeqCst), 1);
        assert!(matches!(result, Err(NistError::NetworkError(err)) if err.is_timeout()));
    }

    #[test]
    fn test_get_data_batch_reports_partial_success() {
        let success_body = http_response(
            "200 OK",
            r#"
            <html>
              <body>
                <h1>Example substance</h1>
                <ol>
                  <li><a href="/cgi/cbook.cgi?ID=example-entry">Example result</a></li>
                </ol>
                <a href="/gas">Gas phase thermochemistry data</a>
                <li>Molecular weight: 16.043</li>
              </body>
            </html>
            "#,
        );
        let (base_url, hits) = spawn_threshold_server(success_body, 4);
        let client = build_timeout_client(250);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let requests = vec![
            ("CH4".to_string(), SearchType::MolarMass, Phase::Gas),
            ("H2O".to_string(), SearchType::MolarMass, Phase::Gas),
        ];

        let report = parser
            .get_data_batch(&requests)
            .expect("batch should partially succeed");

        assert!(hits.load(Ordering::SeqCst) >= 5);
        assert_eq!(report.entries.len(), 2);
        assert_eq!(report.successes(), 1);
        assert_eq!(report.failures(), 1);
        assert!(report.entries[0].result.is_ok());
        assert!(report.entries[1].result.is_err());
        let _ = base_url;
    }

    #[test]
    fn test_get_data_batch_fails_when_every_request_fails() {
        let bad_body = http_response("404 Not Found", "missing");
        let (base_url, hits) = spawn_threshold_server(bad_body, 0);
        let client = build_timeout_client(250);
        let parser =
            NistParser::with_client_and_base_url(client, url::Url::parse(&base_url).unwrap());
        let requests = vec![
            ("CH4".to_string(), SearchType::MolarMass, Phase::Gas),
            ("H2O".to_string(), SearchType::MolarMass, Phase::Gas),
        ];

        let result = parser.get_data_batch(&requests);

        assert!(hits.load(Ordering::SeqCst) >= 1);
        assert!(matches!(
            result,
            Err(NistError::BatchAllFailed { attempts, .. }) if attempts == requests.len()
        ));
        let _ = base_url;
    }

    #[test]
    #[ignore = "live NIST integration test"]
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
    #[ignore = "live NIST integration test"]
    fn test_nonexistent_substance() {
        let parser = NistParser::new();
        let result = parser.get_data("NonexistentSubstance123", SearchType::MolarMass, Phase::Gas);
        assert!(matches!(result, Err(NistError::SubstanceNotFound)));
    }

    #[test]
    #[ignore = "live NIST integration test"]
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
    #[ignore = "live NIST integration test"]
    fn test_simple_substance() {
        let parser = NistParser::new();

        // Test enthalpy data for water
        let result = parser.get_data("O2", SearchType::All, Phase::Gas);
        let dh = result.as_ref().unwrap().clone().dh.unwrap();
        println!("dh: {}", dh);
        assert_eq!(dh, 0.0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_data_structure_initialization() {
        let data = NistInput {
            cp: Some(vec![vec![1.0, 2.0, 3.0]]),
            T: Some(vec![vec![298.0, 1000.0]]),
            dh: Some(42.0),
            ds: Some(100.0),
            molar_mass: Some(18.015),
        };

        assert_eq!(data.cp.as_ref().unwrap()[0][0], 1.0);
        assert_eq!(data.dh.unwrap(), 42.0);
        assert_eq!(data.ds.unwrap(), 100.0);
        assert_eq!(data.molar_mass.unwrap(), 18.015);
    }

    #[test]
    #[ignore = "live NIST integration test"]
    fn test_ch4_gas() {
        let parser = NistParser::new();
        let substance = "CH4";
        let result = parser.get_data(substance, SearchType::All, Phase::Gas);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.cp.is_some());
            assert!(data.dh.is_some());
            assert!(data.ds.is_some());
            assert!(data.molar_mass.is_some());
            println!("CH4 Gas data parsed successfully");
        }
    }

    #[test]
    #[ignore = "live NIST integration test"]
    fn test_nacl_solid() {
        let parser = NistParser::new();
        let substance = "NaCl";
        let result = parser.get_data(substance, SearchType::All, Phase::Solid);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.cp.is_some());
            assert!(data.dh.is_some());
            assert!(data.ds.is_some());
            assert!(data.molar_mass.is_some());
            println!("NaCl Solid data parsed successfully");
        }
    }

    #[test]
    #[ignore = "live NIST integration test"]
    fn test_nacl_liquid() {
        let parser = NistParser::new();
        let substance = "NaCl";
        let result = parser.get_data(substance, SearchType::All, Phase::Liquid);
        assert!(result.is_ok());
        if let Ok(data) = result {
            assert!(data.cp.is_some());
            assert!(data.dh.is_some());
            assert!(data.ds.is_some());
            assert!(data.molar_mass.is_some());
            println!("NaCl Liquid data parsed successfully");
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
