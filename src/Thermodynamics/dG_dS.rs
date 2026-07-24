use crate::Thermodynamics::User_substances::{DataType, Phases, SubsData};
use crate::Thermodynamics::User_substances_error::{SubsDataError, SubsDataResult};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use log::warn;
use std::collections::HashMap;
use std::f64;
const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
const R_sym: Expr = Expr::Const(R);

fn numeric_concentration_at(
    context: &str,
    substance: &str,
    n: &Option<Vec<f64>>,
    index: usize,
) -> Option<f64> {
    let values = n.as_ref()?;
    values.get(index).copied().or_else(|| {
        warn!(
            "Skipping gas correction for '{}' in {} because the concentration vector is too short",
            substance, context
        );
        None
    })
}

fn symbolic_concentration_at(
    context: &str,
    substance: &str,
    n: &Option<Vec<Expr>>,
    index: usize,
) -> Option<Expr> {
    let values = n.as_ref()?;
    values.get(index).cloned().or_else(|| {
        warn!(
            "Skipping symbolic gas correction for '{}' in {} because the concentration vector is too short",
            substance, context
        );
        None
    })
}

fn validate_concentration_vector_length<T>(
    context: &str,
    n: &Option<Vec<T>>,
    expected_len: usize,
) -> SubsDataResult<()> {
    if let Some(values) = n.as_ref() {
        if values.len() < expected_len {
            return Err(SubsDataError::calculation_failed(
                "mixture",
                context,
                format!(
                    "concentration vector length {} is shorter than the number of substances {}",
                    values.len(),
                    expected_len
                ),
            ));
        }
    }
    Ok(())
}

fn with_pure_working_copy<T>(
    subs: &SubsData,
    work: impl FnOnce(&mut SubsData) -> SubsDataResult<T>,
) -> SubsDataResult<T> {
    let mut working_copy = subs.clone();
    work(&mut working_copy)
}

///////////////////////////////////////////////////////////////////dG///////////////////////////////////
impl SubsData {
    /// Attempts to build thermodynamic closures and returns a typed error when
    /// any substance is missing the required calculator state.
    fn calculate_therm_map_of_fun_local(
        &mut self,
    ) -> SubsDataResult<(
        Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>,
        Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>>,
    )> {
        let mut vec_of_cp: Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>> = Vec::new();
        let mut vec_of_dh: Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>> = Vec::new();
        let mut vec_of_ds: Vec<Option<Box<dyn Fn(f64) -> f64 + Send + Sync>>> = Vec::new();
        for substance in &self.substances.clone() {
            let result = self.with_thermo_calculator(substance, |thermo| {
                // Build the closures from the calculator in one place so the
                // lookup state stays canonical.
                thermo.create_closures_Cp_dH_dS()?;
                Ok((
                    Some(thermo.get_C_fun()?),
                    Some(thermo.get_dh_fun()?),
                    Some(thermo.get_ds_fun()?),
                ))
            })?;

            let (cp_fun, dh_fun, ds_fun) = result;
            vec_of_cp.push(cp_fun);
            vec_of_dh.push(dh_fun);
            vec_of_ds.push(ds_fun);
        }
        Ok((vec_of_cp, vec_of_dh, vec_of_ds))
    }

    pub fn calc_dG_for_one_phase(
        &mut self,
        P: f64,

        T: f64,
        n: Option<Vec<f64>>,
        Np: Option<f64>,
    ) -> SubsDataResult<HashMap<String, f64>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length("Gibbs calculation", &n, self.substances.len())?;
        with_pure_working_copy(self, |working| {
            let mut map_to_insert = HashMap::new();
            let Np = Np.unwrap_or(1.0);
            working.calculate_therm_map_of_properties(T)?;
            for (i, substance) in working.substances.iter().enumerate() {
                // Pull numeric values through the explicit readiness API so we do
                // not depend on raw cache shape here.
                let Some(dh) = working
                    .therm_value_state(substance, DataType::dH)
                    .value()
                    .copied()
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dH".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds) = working
                    .therm_value_state(substance, DataType::dS)
                    .value()
                    .copied()
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS".to_string(),
                        substance: substance.clone(),
                    });
                };
                let dG = dh - T * ds;
                // gas_correrction = 0 if w = None or if Phase is not Gas, else gas_correrction = R*T*ln(P/101325) + R*T*ln(w[i])
                // if Phase is Gas or Phase is None
                let gas_correrction: f64 = if let Some(n_i) =
                    numeric_concentration_at("Gibbs calculation", substance, &n, i)
                {
                    let correction: f64 = R * T * f64::ln(P / 101325.0) + R * T * f64::ln(n_i / Np);

                    if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
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
            } // for (i, substance) in self.substances.iter().enumerate()
            Ok(map_to_insert)
        })
    }

    pub fn calculate_Gibbs_sym_one_phase(
        &mut self,
        _T: f64,
        n: Option<Vec<Expr>>,
        Np: Option<Expr>,
    ) -> SubsDataResult<HashMap<String, Expr>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length(
            "symbolic Gibbs calculation",
            &n,
            self.substances.len(),
        )?;
        with_pure_working_copy(self, |working| {
            working.calculate_therm_map_of_sym()?;
            let T = Expr::Var("T".to_owned());
            let P = Expr::Var("P".to_owned());
            let Np = Np.unwrap_or(Expr::Const(1.0));
            let mut map_to_insert = HashMap::new();
            for (i, substance) in working.substances.iter().enumerate() {
                let Some(map_sym) = working.therm_map_of_sym.get(substance).cloned() else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym/dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(dh_sym) = map_sym
                    .get(&DataType::dH_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds_sym) = map_sym
                    .get(&DataType::dS_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let dG_sym = dh_sym - T.clone() * ds_sym;
                let gas_correrction: Expr = if let Some(n_i) =
                    symbolic_concentration_at("symbolic Gibbs calculation", substance, &n, i)
                {
                    let correction: Expr =
                        R_sym * T.clone() * Expr::ln(P.clone() / Expr::Const(101325.0))
                            + R_sym * T.clone() * Expr::ln(n_i / Np.clone());
                    let gas_correrction: Expr = if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
                            _ => Expr::Const(0.0),
                        }
                    } else {
                        correction
                    };
                    gas_correrction
                } else {
                    Expr::Const(0.0)
                };

                let dG_sym = dG_sym.simplify() + gas_correrction.simplify();
                map_to_insert.insert(substance.clone(), dG_sym.clone());
            }
            Ok(map_to_insert)
        })
    }

    pub fn calculate_dG0_sym_one_phase(&mut self) -> SubsDataResult<HashMap<String, Expr>> {
        with_pure_working_copy(self, |working| {
            let T = Expr::Var("T".to_owned());

            working.calculate_therm_map_of_sym()?;
            let mut map_to_insert = HashMap::new();
            for (_, substance) in working.substances.iter().enumerate() {
                let Some(map_sym) = working.therm_map_of_sym.get(substance).cloned() else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym/dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(dh_sym) = map_sym
                    .get(&DataType::dH_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds_sym) = map_sym
                    .get(&DataType::dS_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let dG_sym = dh_sym - T.clone() * ds_sym;
                let dG_sym = dG_sym.simplify();
                map_to_insert.insert(substance.clone(), dG_sym.clone());
            }
            Ok(map_to_insert)
        })
    }
    /// Function for calculating Gibbs free energy of a given mixure of substances ( at given T, P, concentration)
    pub fn calculate_Gibbs_fun_one_phase(
        &mut self,
        P: f64,
        _T: f64,
    ) -> SubsDataResult<
        HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>>,
    > {
        let phases_map = self.map_of_phases.clone();
        with_pure_working_copy(self, |working| {
            let mut map_to_insert = HashMap::new();
            let (_Cp, mut dh_vec, mut ds_vec) = working.calculate_therm_map_of_fun_local()?;
            for (i, substance) in working.substances.iter().enumerate() {
                let Some(dh) = dh_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dH_fun".to_string(),
                    });
                };
                let Some(ds) = ds_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dS_fun".to_string(),
                    });
                };
                let phases = phases_map.clone();

                let gas_correction = {
                    let substance = substance.clone();

                    move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>| -> f64 {
                        if let Some(n_i) =
                            numeric_concentration_at("Gibbs closure calculation", &substance, &n, i)
                        {
                            let Np = Np.unwrap_or(1.0);
                            let correction =
                                R * T * f64::ln(P / 101325.0) + R * T * f64::ln(n_i / Np);
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
                let dG: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static> = dG;
                map_to_insert.insert(substance.clone(), dG);
            }
            Ok(map_to_insert)
        })
    }

    /// Function for calculating Gibbs free energy of a given mixure of substances ( at given T, P, concentration)
    pub fn calculate_dG0_fun_one_phase(
        &mut self,
    ) -> SubsDataResult<HashMap<String, Box<dyn Fn(f64) -> f64 + Send + Sync>>> {
        with_pure_working_copy(self, |working| {
            let mut map_to_insert = HashMap::new();
            let (_Cp, mut dh_vec, mut ds_vec) = working.calculate_therm_map_of_fun_local()?;
            for (i, substance) in working.substances.iter().enumerate() {
                let Some(dh) = dh_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dH_fun".to_string(),
                    });
                };
                let Some(ds) = ds_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dS_fun".to_string(),
                    });
                };

                let dG = Box::new(move |t: f64| dh(t) - t * ds(t));
                let dG: Box<dyn Fn(f64) -> f64 + Send + Sync> = dG;
                map_to_insert.insert(substance.clone(), dG);
            }
            Ok(map_to_insert)
        })
    }
    ////////////////////////////////////////////S - ENTHROPY///////////////////////////////////////////////////////////
    /// Function for calculating enthropy of a given mixure of substances (numerical result at given T, P, concentration)
    pub fn calculate_S_for_one_phase(
        &mut self,
        P: f64,

        T: f64,
        n: Option<Vec<f64>>,
        Np: Option<f64>,
    ) -> SubsDataResult<HashMap<String, f64>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length("entropy calculation", &n, self.substances.len())?;
        with_pure_working_copy(self, |working| {
            let mut map_to_insert = HashMap::new();

            working.calculate_therm_map_of_properties(T)?;
            for (i, substance) in working.substances.iter().enumerate() {
                // The entropy term is read from the explicit state view so the
                // caller cannot confuse a missing cache entry with data.
                let Some(dS) = working
                    .therm_value_state(substance, DataType::dS)
                    .value()
                    .copied()
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS".to_string(),
                        substance: substance.clone(),
                    });
                };

                // gas_correrction = 0 if w = None or if Phase is not Gas, else gas_correrction = R*T*ln(P/101325) + R*T*ln(w[i])
                // if Phase is Gas or Phase is None
                let gas_correrction: f64 = if let Some(n_i) =
                    numeric_concentration_at("entropy calculation", substance, &n, i)
                {
                    let Np = Np.unwrap_or(1.0);
                    let correction: f64 = R * f64::ln(P / 101325.0) + R * f64::ln(n_i / Np);

                    if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
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
            } // for (i, substance) in self.substances.iter().enumerate()
            Ok(map_to_insert)
        })
    }

    pub fn calculate_S_sym_for_one_phase(
        &mut self,
        _T: f64,
        n: Option<Vec<Expr>>,
        Np: Option<Expr>,
    ) -> SubsDataResult<HashMap<String, Expr>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length(
            "symbolic entropy calculation",
            &n,
            self.substances.len(),
        )?;
        with_pure_working_copy(self, |working| {
            working.calculate_therm_map_of_sym()?;
            let P = Expr::Var("P".to_owned());
            let Np = Np.unwrap_or(Expr::Const(1.0));
            let mut map_to_insert = HashMap::new();
            for (i, substance) in working.substances.iter().enumerate() {
                let Some(map_sym) = working.therm_map_of_sym.get(substance).cloned() else {
                    return Err(SubsDataError::MissingData {
                        field: "dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds_sym) = map_sym
                    .get(&DataType::dS_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };

                let gas_correrction: Expr = if let Some(n_i) =
                    symbolic_concentration_at("symbolic entropy calculation", substance, &n, i)
                {
                    let correction: Expr = R_sym * Expr::ln(P.clone() / Expr::Const(101325.0))
                        + R_sym * Expr::ln(n_i / Np.clone());
                    let gas_correrction: Expr = if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
                            _ => Expr::Const(0.0),
                        }
                    } else {
                        correction
                    };
                    gas_correrction
                } else {
                    Expr::Const(0.0)
                };

                let dS_sym = ds_sym.simplify() - gas_correrction.simplify();
                map_to_insert.insert(substance.clone(), dS_sym.clone());
            }
            Ok(map_to_insert)
        })
    }

    pub fn calculate_S_fun_for_one_phase(
        &mut self,
        P: f64,
        _T: f64,
    ) -> SubsDataResult<HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>>
    {
        let phases_map = self.map_of_phases.clone();
        with_pure_working_copy(self, |working| {
            let mut map_to_insert = HashMap::new();
            let (_Cp, _dh_vec, mut ds_vec) = working.calculate_therm_map_of_fun_local()?;
            for (i, substance) in working.substances.iter().enumerate() {
                let Some(ds) = ds_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dS_fun".to_string(),
                    });
                };
                let phases = phases_map.clone();

                let gas_correction = {
                    let substance = substance.clone();

                    move |n: Option<Vec<f64>>, Np: Option<f64>| -> f64 {
                        if let Some(n_i) = numeric_concentration_at(
                            "entropy closure calculation",
                            &substance,
                            &n,
                            i,
                        ) {
                            let Np = Np.unwrap_or(1.0);
                            let correction = R * f64::ln(P / 101325.0) + R * f64::ln(n_i / Np);
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
            Ok(map_to_insert)
        })
    }

    ///////////////////////////////////////////////////////////////////dA - HELMHOLTZ ENERGY///////////////////////////////////
    pub fn calc_dA_for_one_phase(
        &mut self,
        _P: f64,

        T: f64,
        n: Option<Vec<f64>>,
        Np: Option<f64>,
    ) -> SubsDataResult<HashMap<String, f64>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length("Helmholtz calculation", &n, self.substances.len())?;
        with_pure_working_copy(self, |working| {
            working.extract_all_thermal_coeffs(T)?;
            working.calculate_therm_map_of_properties(T)?;
            let mut map_to_insert = HashMap::new();
            let Np = Np.unwrap_or(1.0);

            for (i, substance) in working.substances.iter().enumerate() {
                let Some(dh) = working
                    .therm_value_state(substance, DataType::dH)
                    .value()
                    .copied()
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dH".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds) = working
                    .therm_value_state(substance, DataType::dS)
                    .value()
                    .copied()
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS".to_string(),
                        substance: substance.clone(),
                    });
                };

                // dA = dU - T*dS, where dU = dH - P*V (for ideal gas: P*V = R*T)
                let du = dh - R * T; // For ideal gas
                let dA = du - T * ds;

                let gas_correction: f64 = if let Some(n_i) =
                    numeric_concentration_at("Helmholtz calculation", substance, &n, i)
                {
                    let correction: f64 = R * T * f64::ln(n_i / Np);
                    if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
                            _ => 0.0,
                        }
                    } else {
                        correction
                    }
                } else {
                    0.0
                };

                let dA = dA + gas_correction;
                map_to_insert.insert(substance.clone(), dA);
            }
            Ok(map_to_insert)
        })
    }

    pub fn calculate_Helmholtz_sym_one_phase(
        &mut self,
        T: f64,
        n: Option<Vec<Expr>>,
        Np: Option<Expr>,
    ) -> SubsDataResult<HashMap<String, Expr>> {
        let phases_map = self.map_of_phases.clone();
        validate_concentration_vector_length(
            "symbolic Helmholtz calculation",
            &n,
            self.substances.len(),
        )?;
        with_pure_working_copy(self, |working| {
            working.extract_all_thermal_coeffs(T)?;
            working.calculate_therm_map_of_sym()?;
            let T = Expr::Var("T".to_owned());
            let Np = Np.unwrap_or(Expr::Const(1.0));
            let mut map_to_insert = HashMap::new();

            for (i, substance) in working.substances.iter().enumerate() {
                let Some(map_sym) = working.therm_map_of_sym.get(substance).cloned() else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym/dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(dh_sym) = map_sym
                    .get(&DataType::dH_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dH_sym".to_string(),
                        substance: substance.clone(),
                    });
                };
                let Some(ds_sym) = map_sym
                    .get(&DataType::dS_sym)
                    .and_then(|slot| slot.as_ref())
                    .map(|expr| expr.as_ref().clone())
                else {
                    return Err(SubsDataError::MissingData {
                        field: "dS_sym".to_string(),
                        substance: substance.clone(),
                    });
                };

                // dA = dU - T*dS, where dU = dH - R*T for ideal gas
                let du_sym = dh_sym - R_sym * T.clone();
                let dA_sym = du_sym - T.clone() * ds_sym;

                let gas_correction: Expr = if let Some(n_i) =
                    symbolic_concentration_at("symbolic Helmholtz calculation", substance, &n, i)
                {
                    let correction: Expr = R_sym * T.clone() * Expr::ln(n_i / Np.clone());
                    if let Some(phase) = phases_map.get(substance) {
                        match phase {
                            Some(Phases::Gas) => correction,
                            _ => Expr::Const(0.0),
                        }
                    } else {
                        correction
                    }
                } else {
                    Expr::Const(0.0)
                };

                let dA_sym = dA_sym.simplify() + gas_correction.simplify();
                map_to_insert.insert(substance.clone(), dA_sym.clone());
            }
            Ok(map_to_insert)
        })
    }

    pub fn calculate_Helmholtz_fun_one_phase(
        &mut self,
        _P: f64,
        T: f64,
    ) -> SubsDataResult<HashMap<String, Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64>>>
    {
        let phases_map = self.map_of_phases.clone();
        with_pure_working_copy(self, |working| {
            working.extract_all_thermal_coeffs(T)?;
            let mut map_to_insert = HashMap::new();
            let (_Cp, mut dh_vec, mut ds_vec) = working.calculate_therm_map_of_fun_local()?;

            for (i, substance) in working.substances.iter().enumerate() {
                let Some(dh) = dh_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dH_fun".to_string(),
                    });
                };
                let Some(ds) = ds_vec.get_mut(i).and_then(|slot| slot.take()) else {
                    return Err(SubsDataError::FunctionCreationFailed {
                        substance: substance.clone(),
                        function_type: "dS_fun".to_string(),
                    });
                };
                let phases = phases_map.clone();

                let gas_correction = {
                    let substance = substance.clone();
                    move |T: f64, n: Option<Vec<f64>>, Np: Option<f64>| -> f64 {
                        if let Some(n_i) = numeric_concentration_at(
                            "Helmholtz closure calculation",
                            &substance,
                            &n,
                            i,
                        ) {
                            let Np = Np.unwrap_or(1.0);
                            let correction = R * T * f64::ln(n_i / Np);
                            match phases.get(&substance) {
                                Some(Some(Phases::Gas)) => correction,
                                _ => 0.0,
                            }
                        } else {
                            0.0
                        }
                    }
                };

                let dA = Box::new(move |t: f64, n: Option<Vec<f64>>, Np: Option<f64>| {
                    let du = dh(t) - R * t; // dU = dH - R*T for ideal gas
                    du - t * ds(t) + gas_correction(t, n, Np)
                });
                let dA: Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64> = dA;
                map_to_insert.insert(substance.clone(), dA);
            }
            Ok(map_to_insert)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use RustedSciThe::symbolic::symbolic_engine::Expr;

    #[test]
    fn test_gibbs_paths_return_typed_errors_for_unknown_substances() {
        let mut subs = SubsData::new();
        subs.set_substances(vec!["__definitely_missing__".to_string()]);

        let gibbs_numeric = subs.calc_dG_for_one_phase(101325.0, 298.15, None, None);
        assert!(gibbs_numeric.is_err());

        let gibbs_fun = subs.calculate_Gibbs_fun_one_phase(101325.0, 298.15);
        assert!(gibbs_fun.is_err());

        let gibbs_ref = subs.calculate_dG0_fun_one_phase();
        assert!(gibbs_ref.is_err());

        let entropy_fun = subs.calculate_S_fun_for_one_phase(101325.0, 298.15);
        assert!(entropy_fun.is_err());

        let helmholtz_fun = subs.calculate_Helmholtz_fun_one_phase(101325.0, 298.15);
        assert!(helmholtz_fun.is_err());
    }

    #[test]
    fn test_gibbs_paths_reject_short_concentration_vectors_explicitly() {
        let mut subs = SubsData::new();
        subs.set_substances(vec!["__definitely_missing__".to_string()]);

        let gibbs_numeric = subs.calc_dG_for_one_phase(101325.0, 298.15, Some(vec![]), None);
        assert!(matches!(
            gibbs_numeric,
            Err(SubsDataError::CalculationFailed { .. })
        ));

        let gibbs_symbolic =
            subs.calculate_Gibbs_sym_one_phase(298.15, Some(vec![]), Some(Expr::Const(1.0)));
        assert!(matches!(
            gibbs_symbolic,
            Err(SubsDataError::CalculationFailed { .. })
        ));
    }

    #[test]
    fn test_concentration_helpers_handle_short_vectors_without_panicking() {
        let n_numeric = Some(vec![0.25]);
        assert_eq!(
            numeric_concentration_at("unit-test", "CO2", &n_numeric, 0),
            Some(0.25)
        );
        assert_eq!(
            numeric_concentration_at("unit-test", "CO2", &n_numeric, 1),
            None
        );

        let n_symbolic = Some(vec![Expr::Const(0.25)]);
        assert!(symbolic_concentration_at("unit-test", "CO2", &n_symbolic, 1).is_none());
    }
}
