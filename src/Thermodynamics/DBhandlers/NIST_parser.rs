use RustedSciThe::symbolic::symbolic_engine::Expr;
use prettytable::{Cell, Row, Table};
use reqwest::blocking::Client;
use scraper::{Html, Selector};
use serde::{Deserialize, Serialize};
use thiserror::Error;
use url::Url;
#[allow(non_upper_case_globals, non_snake_case)]
const e2: Expr = Expr::Const(2.0);
#[allow(non_upper_case_globals, non_snake_case)]
const e3: Expr = Expr::Const(3.0);
#[allow(non_upper_case_globals, non_snake_case)]
const e4: Expr = Expr::Const(4.0);
/// HTTP client trait for dependency injection
pub trait HttpClient {
    fn get_text(&self, url: &str) -> Result<String, reqwest::Error>;
}

// Implementation for the real reqwest client
impl HttpClient for Client {
    fn get_text(&self, url: &str) -> Result<String, reqwest::Error> {
        self.get(url).send()?.text()
    }
}
/// error types for the reqwest client
#[derive(Debug, Error)]
#[allow(dead_code)]
pub enum NistError {
    #[error("Network error: {0}")]
    NetworkError(#[from] reqwest::Error),
    #[error("URL parsing error: {0}")]
    UrlError(#[from] url::ParseError),
    #[error("Substance not found")]
    SubstanceNotFound,
    #[error("Invalid data format")]
    InvalidDataFormat,
}
/// struct for the substance data parsed from the NIST web page
#[derive(Debug, Serialize, Deserialize, Clone)]
#[allow(non_snake_case)]
pub struct NistInput {
    pub cp: Option<Vec<Vec<f64>>>,
    pub T: Option<Vec<Vec<f64>>>,
    pub dh: Option<f64>,
    pub coeffs: Option<(f64, f64, f64, f64,f64, f64, f64, f64 )>,
    ds: Option<f64>,
    molar_mass: Option<f64>,
    unit: Option<String>,
    pub unit_multiplier: f64,
}
/// Phase enum: solid, liquid, gas
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum Phase {
    Gas,
    Solid,
    Liquid,
}

#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum SearchType {
    Cp,
    DeltaH,
    DeltaS,
    MolarMass,
    All,
}

impl Phase {
    fn as_str(&self) -> &'static str {
        match self {
            Phase::Gas => "gas",
            Phase::Solid => "solid",
            Phase::Liquid => "liquid",
        }
    }
}

pub struct NistParser<C: HttpClient> {
    client: C,
}

impl NistParser<Client> {
    pub fn new() -> Self {
        Self {
            client: Client::new(),
        }
    }
}

impl<C: HttpClient> NistParser<C> {
    #[allow(dead_code)]
    pub fn with_client(client: C) -> Self {
        Self { client }
    }
    ///////////////////////////////////TRAVELLING THE WEBSITE: http://webbook.nist.gov/cgi/cbook.cgi///////////////////////////////////////////
    pub fn get_data(
        &self,
        substance: &str,
        search_type: SearchType,
        phase: Phase,
    ) -> Result<NistInput, NistError> {
        let url = self.construct_url(substance)?;
        println!("\n \n URL found: {}", url);
        let html = self.fetch_page(&url)?;
        if !self.check_substance_exists(&html) {
            return Err(NistError::SubstanceNotFound);
        }

        let url_of_substance = self.get_url_of_substance(&html, &url)?;
        println!(
            "\n \n URL of substance {} found: {}",
            substance, url_of_substance
        );
        let html_of_substance = self.fetch_page(&url_of_substance)?;

        let final_url = self.get_final_url(&html_of_substance, &url_of_substance, phase)?;
        let html_of_phase = self.fetch_page(&final_url)?;
        println!("\n \n Final URL found: {}", final_url);
        //   println!("\n \n HTML of phase: {}", html_of_phase);

        let data = self.parse_data(&final_url, search_type, phase)?;
        Ok(data)
    }

    fn construct_url(&self, substance: &str) -> Result<Url, NistError> {
        let substance = substance.replace(' ', "");

        // Try to determine if it's a CAS number (contains '-')
        if substance.contains('-') {
            Ok(Url::parse(&format!(
                "https://webbook.nist.gov/cgi/cbook.cgi?ID={}&Units=SI",
                substance
            ))?)
        } else if substance.chars().any(|c| c.is_ascii_digit()) {
            // If contains numbers, assume it's a chemical formula
            Ok(Url::parse(&format!(
                "https://webbook.nist.gov/cgi/cbook.cgi?Formula={}&NoIon=on&Units=SI",
                substance
            ))?)
        } else {
            // Otherwise, assume it's a name
            Ok(Url::parse(&format!(
                "https://webbook.nist.gov/cgi/cbook.cgi?Name={}&Units=SI",
                substance
            ))?)
        }
    }

    fn fetch_page(&self, url: &Url) -> Result<String, NistError> {
        Ok(self.client.get_text(url.as_str())?)
    }

    fn check_substance_exists(&self, html: &str) -> bool {
        let document = Html::parse_document(html);
        let selector = Selector::parse("h1").unwrap();

        for element in document.select(&selector) {
            let text = element.text().collect::<String>();
            if text.contains("Not Found") {
                return false;
            }
        }
        true
    }

    fn get_url_of_substance(&self, html: &str, original_url: &Url) -> Result<Url, NistError> {
        let document = Html::parse_document(html);

        // Check if we're on a search results page
        if let Ok(selector) = Selector::parse("ol li a") {
            if let Some(first_result) = document.select(&selector).next() {
                if let Some(href) = first_result.value().attr("href") {
                    return Ok(Url::parse(&format!("https://webbook.nist.gov{}", href))?);
                }
            }
        }

        // If not on search results page, use original URL
        Ok(original_url.clone())
    }

    fn get_final_url(
        &self,
        html: &str,
        url_of_substance: &Url,
        phase: Phase,
    ) -> Result<Url, NistError> {
        let document = Html::parse_document(html);
        let link_text = match phase {
            Phase::Gas => "Gas phase thermochemistry data",
            Phase::Solid => "Condensed phase thermochemistry data",
            Phase::Liquid => "Condensed phase thermochemistry data",
        };

        let selector = Selector::parse("a").unwrap();
        for element in document.select(&selector) {
            if element.text().collect::<String>().contains(link_text) {
                if let Some(href) = element.value().attr("href") {
                    return url_of_substance
                        .join(href)
                        .map_err(|e| NistError::UrlError(e));
                }
            }
        }

        // If we couldn't find the link, return the original URL
        Ok(url_of_substance.clone())
    }

    ////////////////////////////////PARSING DATA//////////////////////////////////////////////////////////////
    fn parse_data(
        &self,
        url: &Url,
        search_type: SearchType,
        phase: Phase,
    ) -> Result<NistInput, NistError> {
        print!("\n \n Fetching data...");
        let html = self.fetch_page(url)?;
        let document = Html::parse_document(&html);
        //  println!("\n \n document: {:?}", document);

        let mut data = NistInput {
            cp: None,
            T: None,
            dh: None,
            ds: None,
            coeffs: None,
            molar_mass: None,
            unit: None,
            unit_multiplier: 1.0,
        };

        match search_type {
            SearchType::Cp => {
                (data.cp, data.T) = self.extract_cp(&document, phase)?;
            }
            SearchType::DeltaH | SearchType::DeltaS => {
                let (dh, ds) = self.extract_thermodynamic_data(&document, phase)?;
                data.dh = dh;
                data.ds = ds;
            }
            SearchType::MolarMass => {
                data.molar_mass = self.extract_molar_mass(&document)?;
            }
            SearchType::All => {
                (data.cp, data.T) = self.extract_cp(&document, phase)?;
                let (dh, ds) = self.extract_thermodynamic_data(&document, phase)?;
                data.dh = dh;
                data.ds = ds;
                data.molar_mass = self.extract_molar_mass(&document)?;
            }
        }
        println!("\n \n Data found: {:?} \n \n", data);
        Ok(data)
    }

    fn extract_cp(
        &self,
        document: &Html,
        phase: Phase,
    ) -> Result<(Option<Vec<Vec<f64>>>, Option<Vec<Vec<f64>>>), NistError> {
        let table_selector = match phase {
            //matches the phase parameter to determine the appropriate CSS selector for the table.
            Phase::Gas => {
                Selector::parse("table[aria-label='Gas Phase Heat Capacity (Shomate Equation)']")
                    .unwrap()
            }
            Phase::Solid => {
                Selector::parse("table[aria-label='Solid Phase Heat Capacity (Shomate Equation)']")
                    .unwrap()
            }
            Phase::Liquid => {
                Selector::parse("table[aria-label='Liquid Phase Heat Capacity (Shomate Equation)']")
                    .unwrap()
            }
        }; //Selector::parse(...).unwrap() parses the CSS selector string and unwraps the result, assuming it is valid.

        if let Some(table) = document.select(&table_selector).next() {
            //attempt to select the first table element that
            //matches the table_selector from the document. document.select(&table_selector) returns an iterator over matching elements.

            println!("\n \n found table: {:?} \n \n", table);
            #[allow(non_snake_case)]
            let mut headers_T: Vec<Vec<f64>> = Vec::new();
            // Parse headers
            if let Some(header_row) = table.select(&Selector::parse("tr").unwrap()).next() {
                headers_T = header_row
                    .select(&Selector::parse("td").unwrap())
                    .filter_map(|cell| {
                        let text = cell.text().collect::<String>();

                        let temps: Vec<f64> = text
                            .split_whitespace() // Split the string by whitespace
                            .filter_map(|s| s.parse::<f64>().ok()) // Parse each part to f32, ignoring any errors
                            .collect(); // Collect the results into a Vec<f32>

                        Some(temps)
                    })
                    .collect();
            }
            //    println!("\n \n headers: {:?} \n \n", headers_T);
            let mut coefficients: Vec<Vec<f64>> = vec![Vec::with_capacity(9); headers_T.len()];
            for row in table.select(&Selector::parse("tr").unwrap()).skip(1) {
                //skip the first row, which contains the header
                let cells: Vec<f64> = row
                    .select(&Selector::parse("td").unwrap())
                    .filter_map(|cell| {
                        cell.text()
                            .collect::<String>()
                            .split_whitespace()
                            .next()
                            .and_then(|s| s.parse::<f64>().ok())
                    })
                    .collect();

                if !cells.is_empty() {
                    for (i, cell) in cells.iter().enumerate() {
                        coefficients[i].push(*cell);
                    }
                }
            }

            if !coefficients.is_empty() && !headers_T.is_empty() {
                print!("Cp(T) parsed");
                return Ok((Some(coefficients), Some(headers_T)));
            }
        }

        Ok((None, None))
    }

    fn extract_thermodynamic_data(
        &self,
        document: &Html,
        phase: Phase,
    ) -> Result<(Option<f64>, Option<f64>), NistError> {
        let table_selector = Selector::parse("table").unwrap();

        if let Some(table) = document.select(&table_selector).next() {
            println!("\n \n table: {:?} \n \n", table);
            let mut dh = None;
            let mut ds = None;

            for row in table.select(&Selector::parse("tr").unwrap()) {
                let cells: Vec<String> = row
                    .select(&Selector::parse("td").unwrap())
                    .map(|cell| cell.text().collect::<String>())
                    .collect();

                if cells.len() >= 2 {
                    //this is a numerical value of dH or dS
                    let value_str = cells[1].trim();
                    if let Some(value) = value_str.split('±').next() {
                        if cells[0].contains("H°")
                            && cells[0].contains("f")
                            && cells[0].contains(phase.as_str())
                        {
                            //    println!("\n \n value: {:?} \n \n", value);
                            if let Ok(val) = value.trim().parse::<f64>() {
                                dh = Some(val);
                            }
                        } else if cells[0].starts_with("S°") && cells[0].contains(phase.as_str()) {
                            if let Ok(val) = value.trim().parse::<f64>() {
                                ds = Some(val);
                            }
                        }
                    }
                }
            }

            Ok((dh, ds))
        } else {
            Ok((None, None))
        }
    }

    fn extract_molar_mass(&self, document: &Html) -> Result<Option<f64>, NistError> {
        let selector = Selector::parse("li").unwrap();

        for element in document.select(&selector) {
            let text = element.text().collect::<String>();
            if text.contains("Molecular weight") {
                if let Some(value) = text.split(':').nth(1) {
                    if let Ok(mass) = value.trim().parse::<f64>() {
                        return Ok(Some(mass));
                    }
                }
            }
        }

        Ok(None)
    }
}

impl NistInput {
    pub fn new() -> NistInput {
        NistInput {
            cp: None,
            T: None,
            dh: None,
            ds: None,
            coeffs: None,
            molar_mass: None,
            unit: None,
            unit_multiplier: 1.0,
        }
    }

    pub fn set_unit(&mut self, unit: &str) {
        match unit {
            "J" => {
                self.unit = Some("J".to_string());
                self.unit_multiplier = 1.0;
            }
            "cal" => {
                self.unit = Some("cal".to_string());
                self.unit_multiplier = 1.0 / 4.184;
            }
            _ => panic!("Invalid unit: {}", unit),
        }
    }
    /// Prints a pretty table of the parsed data to the console. The table shows the heat capacity coefficients for every temperature range, and the molar mass, standard enthalpy of formation, and standard entropy of the substance.

    pub fn pretty_print(&self) {
        let temps = self.T.clone().expect("no temperature parsed");
        let coeffs = self.cp.clone().expect("no coefficients parsed");
        let mut table = Table::new();
        let mut header_row = vec![Cell::new("Coefficients")];
        for T in temps.clone() {
            let Tstr = format!("{:} - {:}   ", T[0], T[1]);
            header_row.push(Cell::new(&Tstr));
        }
        table.add_row(Row::new(header_row));

        let coeffs_names = vec!["A", "B", "C", "D", "E", "F", "G", "H"];
        for (i, coeff_name) in coeffs_names.iter().enumerate() {
            let mut row = vec![Cell::new(coeff_name)];

            for coeff_for_every_T in coeffs.iter() {
                let coeff = coeff_for_every_T[i].to_string();
                row.push(Cell::new(&format!("{:}", coeff)));
            }
            table.add_row(Row::new(row));
        }

        table.printstd();

        let mut table = Table::new();
        let M = self.molar_mass.clone().expect("no molar mass parsed");
        let dH = self.dh.clone().expect("no dH parsed");
        let dS = self.ds.clone().expect("no dS parsed");
        let header_row = vec![Cell::new("Molar mass"), Cell::new("dH"), Cell::new("dS")];
        table.add_row(Row::new(header_row));
        let row = vec![
            Cell::new(&M.to_string()),
            Cell::new(&dH.to_string()),
            Cell::new(&dS.to_string()),
        ];
        table.add_row(Row::new(row));
        table.printstd();
    }

    pub fn extract_coefficients(&mut self, T: f64) -> Result<(f64, f64, f64, f64,f64, f64, f64, f64 ), std::io::Error> {
    
        for (i, T_pairs) in self.T.clone().unwrap().iter().enumerate() {
            if T >= T_pairs[0] && T <= T_pairs[1] {
          
                let coeffs = self.cp.clone().unwrap()[i].clone();
                
                let (a, b, c, d, e, f, g, h) = (
                    coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5], coeffs[6],
                    coeffs[7],
                );
                self.coeffs = Some((a, b, c, d, e, f, g, h));
                return Ok((a, b, c, d, e, f, g, h));
            }
        }

        Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "No temperature range found for the given temperature",
        ))
        
    }
    pub fn caclc_cp_dh_ds(&self, T: f64) -> Result<(f64, f64, f64), std::io::Error> {
        let um = self.unit_multiplier;
            if let Some(coeffs )  = self.coeffs.clone(){
                let T = T / 1000.0;
                let (a, b, c, d, e, f, g, h) = (
                    coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6,
                    coeffs.7,
                );
                let Cp = um * calculate_cp(T, a, b, c, d, e);
                let dh0 = self.dh.clone().unwrap();
                //  let ds0 = self.ds.clone().unwrap();
                let dh = um * (calculate_dh(T, a, b, c, d, e, f, g, h) + dh0);
                let ds = um * calculate_s(T, a, b, c, d, e, f, g, h);
                println!(
                    "\n \n Cp: {:?} \n \n dh: {:?} \n \n ds: {:?} \n \n",
                    Cp, dh, ds
                );
                return Ok((Cp, dh, ds));
    
            }
        Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "No temperature range found for the given temperature",
        ))
    }

    pub fn create_sym_cp_dh_ds(&self) -> Result<(Expr, Expr, Expr), std::io::Error> {
        let um = Expr::Const(self.unit_multiplier);
        if let Some(coeffs )  = self.coeffs.clone(){
            let (a, b, c, d, e, f, g, h) = (
                coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6,
                coeffs.7,);
                let Cp = um.clone() * calculate_cp_sym(a, b, c, d, e);
                let dh0 = self.dh.clone().unwrap();
                let dh = um.clone() * (calculate_dh_sym(a, b, c, d, e, f, g, h) + Expr::Const(dh0));
                let ds = um.clone() * calculate_s_sym(a, b, c, d, e, f, g, h);
                println!(
                    "\n \n Cp: {:?} \n \n dh: {:?} \n \n ds: {:?} \n \n",
                    Cp, dh, ds
                );
                return Ok((Cp.symplify(), dh.symplify(), ds.symplify()));
            }

        Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "No temperature range found for the given temperature",
        ))
    }

    pub fn create_closure_cp_dh_ds(
        &self,
    
    ) -> Result<
        (
            Box<dyn Fn(f64) -> f64>,
            Box<dyn Fn(f64) -> f64>,
            Box<dyn Fn(f64) -> f64>,
        ),
        std::io::Error,
    > {
       
            let um = self.unit_multiplier;
            if let Some(coeffs )  = self.coeffs.clone(){
                let (a, b, c, d, e, f, g, h) = (
                    coeffs.0, coeffs.1, coeffs.2, coeffs.3, coeffs.4, coeffs.5, coeffs.6,
                    coeffs.7,
                );
                let Cp = Box::new(move |t| um * calculate_cp(t / 1000.0, a, b, c, d, e));
                let dh0 = self.dh.clone().unwrap();
                let dh = Box::new(move |t| {
                    um * (calculate_dh(t / 1000.0, a, b, c, d, e, f, g, h) + dh0)
                });
                let ds = Box::new(move |t| um * calculate_s(t / 1000.0, a, b, c, d, e, f, g, h));
                return Ok((Cp, dh, ds));
            }
        

        Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "No temperature range found for the given temperature",
        ))
    }
}
////////////////////////////////////////NIST FORMAT FUNCTIONS//////////////////////////////////////////////////////////////////////////
pub fn calculate_cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    a + b * t + c * t.powi(2) + d * t.powi(3) + e / t.powi(2)
}

pub fn calculate_dh(
    t: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    f: f64,
    _g: f64,
    h: f64,
) -> f64 {
    a * t + (b * t.powi(2)) / 2.0 + (c * t.powi(3)) / 3.0 + (d * t.powi(4)) / 4.0 - e / t + f - h
}

pub fn calculate_s(
    t: f64,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    _f: f64,
    g: f64,
    _h: f64,
) -> f64 {
    a * t.ln() + b * t + (c * t.powi(2)) / 2.0 + (d * t.powi(3)) / 3.0 - e / (2.0 * t.powi(2)) + g
}

fn calculate_cp_sym(a: f64, b: f64, c: f64, d: f64, e: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a)
        + Expr::Const(b) * t.clone()
        + Expr::Const(c) * t.clone().pow(e2)
        + Expr::Const(d) * t.clone().pow(e3)
        + Expr::Const(e) / t.pow(e2)
}

fn calculate_dh_sym(a: f64, b: f64, c: f64, d: f64, e: f64, f: f64, _g: f64, h: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a) * t.clone()
        + (Expr::Const(b) * t.clone().pow(e2)) / e2
        + (Expr::Const(c) * t.clone().pow(e3)) / e3
        + (Expr::Const(d) * t.clone().pow(e4)) / e4
        - Expr::Const(e) / t.clone()
        + Expr::Const(f)
        - Expr::Const(h)
}

fn calculate_s_sym(a: f64, b: f64, c: f64, d: f64, e: f64, _f: f64, g: f64, _h: f64) -> Expr {
    let T = Expr::Var("T".to_owned());
    let t = T / Expr::Const(1000.0);
    Expr::Const(a) * t.clone().ln()
        + Expr::Const(b) * t.clone()
        + (Expr::Const(c) * t.clone().pow(e2)) / e2
        + (Expr::Const(d) * t.clone().pow(e3)) / e3
        - Expr::Const(e) / (e2 * t.pow(e2))
        + Expr::Const(g)
}
#[cfg(test)]
mod tests {
    use super::*;
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

/*
TODO!
Added a HttpClient trait for dependency injection
Modified NistParser to be generic over the HTTP client type
Added comprehensive tests including:
Unit tests for URL construction
Tests for phase conversion
Tests for data structure initialization
Mock tests for HTTP requests
Error handling tests
The test suite now covers:
URL Construction:
CAS number URLs
Chemical formula URLs
Chemical name URLs
Space handling in names
Phase Conversion:
Gas phase string conversion
Solid phase string conversion
Real Substance Fetching:
Methane data fetch
Water data fetch
Thermodynamic data fetch
Mock Tests:
Successful substance fetch with mocked response
Not found error handling
Network error handling
Data Structure:
Initialization and access of  NistInput fields
*/
