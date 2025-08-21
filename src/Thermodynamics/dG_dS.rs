use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::User_substances::{
    CalculatorType, DataType, Phases, SearchResult, SubsData, WhatIsFound,
};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;
use std::f64;
const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
const R_sym: Expr = Expr::Const(R);

///////////////////////////////////////////////////////////////////dG///////////////////////////////////
pub fn calc_dG_for_one_phase(
    P: f64,
    subdata: &mut SubsData,
    T: f64,
    n: Option<Vec<f64>>,
    Np: Option<f64>,
) -> HashMap<String, f64> {
    let phases_map = subdata.map_of_phases.clone();

    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let _ = sd.calculate_therm_map_of_properties(T);
    let mut map_to_insert = HashMap::new();
    let Np = Np.unwrap_or(1.0);
    for (i, substance) in sd.substances.iter().enumerate() {
        let map_property_values = sd.therm_map_of_properties_values.get(substance).unwrap();

        let dh = map_property_values.get(&DataType::dH).unwrap().unwrap();
        let ds = map_property_values.get(&DataType::dS).unwrap().unwrap();
        let dG = dh - T * ds;
        println!(
            "substance:{}: dh {}, ds {}, dG {} \n",
            substance, dh, ds, dG
        );
        // gas_correrction = 0 if w = None or if Phase is not Gas, else gas_correrction = R*T*ln(P/101325) + R*T*ln(w[i])
        // if Phase is Gas or Phase is None
        let gas_correrction: f64 = if let Some(ref n) = n {
            let correction: f64 = R * T * f64::ln(P / 101325.0) + R * T * f64::ln(n[i] / Np);

            if let Some(phase) = phases_map.get(substance).unwrap() {
                match phase {
                    Phases::Gas => correction,
                    _ => 0.0,
                }
            } else {
                correction
            }
        } else {
            0.0
        };
        let dG = dG + gas_correrction;
        map_to_insert.insert(substance.clone(), dG);
    } // for (i, substance) in sd.substances.iter().enumerate()
    map_to_insert
}

pub fn calculate_Gibbs_sym_one_phase(
    subdata: &mut SubsData,
    T: f64,
    n: Option<Vec<Expr>>,
    Np: Option<Expr>,
) -> HashMap<String, Expr> {
    let phases_map = subdata.map_of_phases.clone();
    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let _ = sd.calculate_therm_map_of_sym();
    let T = Expr::Var("T".to_owned());
    let P = Expr::Var("P".to_owned());
    let Np = Np.unwrap_or(Expr::Const(1.0));
    let mut map_to_insert = HashMap::new();
    for (i, substance) in sd.substances.iter().enumerate() {
        let map_sym = sd
            .clone()
            .therm_map_of_sym
            .get(&substance.clone())
            .unwrap()
            .clone();
        let dh_sym = *map_sym.get(&DataType::dH_sym).unwrap().clone().unwrap();
        let ds_sym = *map_sym.get(&DataType::dS_sym).unwrap().clone().unwrap();
        let dG_sym = dh_sym - T.clone() * ds_sym;
        let gas_correrction: Expr = if let Some(ref n) = n {
            let correction: Expr = R_sym * T.clone() * Expr::ln(P.clone() / Expr::Const(101325.0))
                + R_sym * T.clone() * Expr::ln(n[i].clone() / Np.clone());
            let gas_correrction: Expr = if let Some(phase) = phases_map.get(substance).unwrap() {
                match phase {
                    Phases::Gas => correction,
                    _ => Expr::Const(0.0),
                }
            } else {
                correction
            };
            gas_correrction
        } else {
            Expr::Const(0.0)
        };

        let dG_sym = dG_sym.symplify() + gas_correrction.symplify();
        map_to_insert.insert(substance.clone(), dG_sym.clone());
    }
    map_to_insert
}

/// Function for calculating Gibbs free energy of a given mixure of substances ( at given T, P, concentration)
pub fn calculate_Gibbs_fun_one_phase(
    subdata: &mut SubsData,
    P: f64,
    T: f64,
) -> HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>> {
    let phases_map = subdata.map_of_phases.clone();
    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let mut map_to_insert = HashMap::new();
    let (_Cp, mut dh_vec, mut ds_vec) = calculate_therm_map_of_fun_local(sd).unwrap();
    for (i, substance) in sd.substances.iter().enumerate() {
        let dh = dh_vec[i].take().unwrap();
        let ds = ds_vec[i].take().unwrap();
        let phases = phases_map.clone();

        let gas_correction = {
            let substance = substance.clone();

            move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>| -> f64 {
                if let Some(n) = n {
                    let Np = Np.unwrap_or(1.0);
                    let correction = R * T * f64::ln(P / 101325.0) + R * T * f64::ln(n[i] / Np);
                    match phases.get(&substance) {
                        Some(Some(Phases::Gas)) => correction,
                        _ => 0.0,
                    }
                } else {
                    0.0
                }
            }
        };
        let dG = Box::new(move |t: f64, n: Option<Vec<f64>>, Np: Option<f64>| {
            dh(t) - t * ds(t) + gas_correction(t, n, Np)
        });
        let dG: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64> = dG;
        map_to_insert.insert(substance.clone(), dG);
    }
    map_to_insert
}

/// Creates closures for the thermodynamic properties of the substances
pub fn calculate_therm_map_of_fun_local(
    sd: &mut SubsData,
) -> Result<
    (
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>>,
    ),
    String,
> {
    let mut vec_of_cp: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    let mut vec_of_dh: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    let mut vec_of_ds: Vec<Option<Box<dyn Fn(f64) -> f64 + 'static>>> = Vec::new();
    for substance in &sd.substances.clone() {
        match sd
            .search_results
            .get_mut(substance)
            .unwrap()
            .get_mut(&WhatIsFound::Thermo)
            .unwrap()
        {
            Some(SearchResult {
                calculator: Some(CalculatorType::Thermo(thermo)),
                ..
            }) => {
                let mut thermo = thermo.clone();

                // Create closures for thermodynamic functions
                if let Err(e) = thermo.create_closures_Cp_dH_dS() {
                    println!(
                        "Warning: Failed to create closures for {}: {}",
                        substance, e
                    );
                    continue;
                }

                match thermo.get_C_fun() {
                    Ok(cp_fun) => {
                        vec_of_cp.push(Some(cp_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get Cp function for {}: {}",
                            substance, e
                        );
                    }
                };

                match thermo.get_dh_fun() {
                    Ok(dh_fun) => {
                        vec_of_dh.push(Some(dh_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get dH function for {}: {}",
                            substance, e
                        );
                    }
                };

                match thermo.get_ds_fun() {
                    Ok(ds_fun) => {
                        vec_of_ds.push(Some(ds_fun));
                    }
                    Err(e) => {
                        println!(
                            "Warning: Failed to get dS function for {}: {}",
                            substance, e
                        );
                    }
                };
            }
            _ => {
                continue;
            }
        }
    }
    Ok((vec_of_cp, vec_of_dh, vec_of_ds))
}
////////////////////////////////////////////S - ENTHROPY///////////////////////////////////////////////////////////
/// Function for calculating enthropy of a given mixure of substances (numerical result at given T, P, concentration)
pub fn calculate_S_for_one_phase(
    P: f64,
    subdata: &mut SubsData,
    T: f64,
    n: Option<Vec<f64>>,
    Np: Option<f64>,
) -> HashMap<String, f64> {
    let phases_map = subdata.map_of_phases.clone();

    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let _ = sd.calculate_therm_map_of_properties(T);
    let mut map_to_insert = HashMap::new();

    for (i, substance) in sd.substances.iter().enumerate() {
        let map_property_values = sd.therm_map_of_properties_values.get(substance).unwrap();

        let dS = map_property_values.get(&DataType::dS).unwrap().unwrap();

        // gas_correrction = 0 if w = None or if Phase is not Gas, else gas_correrction = R*T*ln(P/101325) + R*T*ln(w[i])
        // if Phase is Gas or Phase is None
        let gas_correrction: f64 = if let Some(ref n) = n {
            let Np = Np.unwrap_or(1.0);
            let correction: f64 = R * f64::ln(P / 101325.0) + R * f64::ln(n[i] / Np);

            if let Some(phase) = phases_map.get(substance).unwrap() {
                match phase {
                    Phases::Gas => correction,
                    _ => 0.0,
                }
            } else {
                correction
            }
        } else {
            0.0
        };
        let dS = dS - gas_correrction;
        map_to_insert.insert(substance.clone(), dS);
    } // for (i, substance) in sd.substances.iter().enumerate()
    map_to_insert
}

pub fn calculate_S_sym_for_one_phase(
    subdata: &mut SubsData,
    T: f64,
    n: Option<Vec<Expr>>,
    Np: Option<Expr>,
) -> HashMap<String, Expr> {
    let phases_map = subdata.map_of_phases.clone();
    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let _ = sd.calculate_therm_map_of_sym();

    let P = Expr::Var("P".to_owned());
    let Np = Np.unwrap_or(Expr::Const(1.0));
    let mut map_to_insert = HashMap::new();
    for (i, substance) in sd.substances.iter().enumerate() {
        let map_sym = sd
            .clone()
            .therm_map_of_sym
            .get(&substance.clone())
            .unwrap()
            .clone();

        let ds_sym = *map_sym.get(&DataType::dS_sym).unwrap().clone().unwrap();

        let gas_correrction: Expr = if let Some(ref n) = n {
            let correction: Expr = R_sym * Expr::ln(P.clone() / Expr::Const(101325.0))
                + R_sym * Expr::ln(n[i].clone() / Np.clone());
            let gas_correrction: Expr = if let Some(phase) = phases_map.get(substance).unwrap() {
                match phase {
                    Phases::Gas => correction,
                    _ => Expr::Const(0.0),
                }
            } else {
                correction
            };
            gas_correrction
        } else {
            Expr::Const(0.0)
        };

        let dS_sym = ds_sym.symplify() - gas_correrction.symplify();
        map_to_insert.insert(substance.clone(), dS_sym.clone());
    }
    map_to_insert
}

pub fn calculate_S_fun_for_one_phase(
    subdata: &mut SubsData,
    P: f64,
    T: f64,
) -> HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>> {
    let phases_map = subdata.map_of_phases.clone();
    let sd = subdata;
    let _ = sd.extract_all_thermal_coeffs(T);
    let mut map_to_insert = HashMap::new();
    let (_Cp, _dh_vec, mut ds_vec) = calculate_therm_map_of_fun_local(sd).unwrap();
    for (i, substance) in sd.substances.iter().enumerate() {
        let ds = ds_vec[i].take().unwrap();
        let phases = phases_map.clone();

        let gas_correction = {
            let substance = substance.clone();

            move |n: Option<Vec<f64>>, Np: Option<f64>| -> f64 {
                if let Some(n) = n {
                    let Np = Np.unwrap_or(1.0);
                    let correction = R * f64::ln(P / 101325.0) + R * f64::ln(n[i] / Np);
                    match phases.get(&substance) {
                        Some(Some(Phases::Gas)) => correction,
                        _ => 0.0,
                    }
                } else {
                    0.0
                }
            }
        };
        let dS = Box::new(move |t: f64, n: Option<Vec<f64>>, Np: Option<f64>| {
            ds(t) - gas_correction(n, Np)
        });
        let dS: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64> = dS;
        map_to_insert.insert(substance.clone(), dS);
    }
    map_to_insert
}
