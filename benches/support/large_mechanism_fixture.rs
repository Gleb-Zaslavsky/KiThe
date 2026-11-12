// Shared physical large-mechanism fixture for correctness stories and Criterion.

use std::collections::{BTreeMap, HashMap};

use nalgebra::DMatrix;
use root::Kinetics::molmass::calculate_molar_mass;
use root::SimpleReactorsBVP::SimpleReactorBVP::{FastElemReact, R_G, SimpleReactorTask};
use root::SimpleReactorsBVP::reactor_BVP_utils::ScalingConfig;
use root::Thermodynamics::User_substances::{LibraryPriority, SubsData, WhatIsFound};

pub const REACTION_IDS: &[&str] = &["296", "297", "299", "300", "301"];
pub const THERMO_TEMPERATURE_K: f64 = 1200.0;
pub const REACTOR_LENGTH_M: f64 = 0.01;
const TEMPERATURE_SCALE_K: f64 = 300.0;
const MASS_FLUX_KG_PER_M2_S: f64 = 0.2;
const GAS_DENSITY_KG_PER_M3: f64 = 0.03;
const THERMAL_CONDUCTIVITY_W_PER_M_K: f64 = 0.05;
const SPECIES_DIFFUSIVITY_M2_PER_S: f64 = 2.0e-5;
const INLET_MOLE_FRACTIONS: [(&str, f64); 8] = [
    ("CH3CHO", 0.012),
    ("O2", 0.18),
    ("N2", 0.78498),
    ("H2O", 0.02),
    ("OH", 0.001),
    ("H", 1.0e-5),
    ("O", 1.0e-5),
    ("HO2", 0.002),
];

#[derive(Debug)]
pub struct ReactionFixture {
    pub equation: String,
    pub a_si: f64,
    pub temperature_exponent: f64,
    pub activation_energy_j_per_mol: f64,
}

fn parse_side(side: &str) -> Vec<(String, f64)> {
    side.split('+')
        .map(str::trim)
        .filter(|term| !term.is_empty())
        .map(|term| {
            let mut words = term.split_whitespace();
            let first = words.next().expect("nonempty reaction term");
            match first.parse::<f64>() {
                Ok(coefficient) => (
                    words
                        .next()
                        .expect("species after stoichiometric coefficient")
                        .to_string(),
                    coefficient,
                ),
                Err(_) => (term.to_string(), 1.0),
            }
        })
        .collect()
}

pub fn parse_equation(equation: &str) -> (Vec<(String, f64)>, Vec<(String, f64)>) {
    let (reactants, products) = equation
        .split_once("=>")
        .expect("fixture contains only one-way elementary reactions");
    (parse_side(reactants), parse_side(products))
}

pub fn load_reactions() -> (Vec<ReactionFixture>, Vec<String>) {
    let root: serde_json::Value = serde_json::from_str(include_str!("../../Reactbase.json"))
        .expect("bundled reaction database is valid JSON");
    let cantera = root
        .get("Cantera")
        .and_then(serde_json::Value::as_object)
        .expect("bundled Cantera mechanism exists");
    let mut species = BTreeMap::<String, ()>::new();
    let mut reactions = Vec::with_capacity(REACTION_IDS.len());

    for id in REACTION_IDS {
        let row = cantera
            .get(*id)
            .unwrap_or_else(|| panic!("missing bundled Cantera reaction {id}"));
        let equation = row["eq"].as_str().expect("reaction equation").to_string();
        assert!(equation.contains("=>") && !equation.contains("<=>"));
        assert_eq!(row["type"].as_str(), Some("elem"));
        let arrhenius = row["Arrenius"].as_array().expect("Arrhenius triple");
        assert_eq!(arrhenius.len(), 3);
        let source_a = arrhenius[0].as_f64().expect("finite A");
        let temperature_exponent = arrhenius[1].as_f64().expect("finite n");
        let activation_energy_j_per_mol = arrhenius[2].as_f64().expect("finite E");
        let (reactants, products) = parse_equation(&equation);
        let order: f64 = reactants.iter().map(|(_, coefficient)| coefficient).sum();
        let a_si = source_a * 1.0e6_f64.powf(1.0 - order);
        for (name, _) in reactants.iter().chain(&products) {
            species.insert(name.clone(), ());
        }
        reactions.push(ReactionFixture {
            equation,
            a_si,
            temperature_exponent,
            activation_energy_j_per_mol,
        });
    }

    species.insert("N2".to_string(), ());
    (reactions, species.into_keys().collect())
}

pub fn resolved_thermochemistry(
    species: &[String],
) -> (
    HashMap<String, (f64, f64)>,
    HashMap<String, f64>,
    HashMap<String, String>,
) {
    let mut subs = SubsData::new();
    subs.set_substances(species.to_vec());
    subs.set_multiple_library_priorities(
        vec![
            "NASA_gas".to_string(),
            "NASA".to_string(),
            "NASA7".to_string(),
            "nuig_thermo".to_string(),
        ],
        LibraryPriority::Priority,
    );
    subs.search_substances()
        .expect("local thermochemistry lookup should complete");
    let missing = subs.get_not_found_substances();
    assert!(
        missing.is_empty(),
        "missing bundled thermo records: {missing:?}"
    );

    let mut thermo = HashMap::new();
    let mut molar_masses = HashMap::new();
    let mut provenance = HashMap::new();
    for name in species {
        subs.extract_thermal_coeffs(name, THERMO_TEMPERATURE_K)
            .unwrap_or_else(|error| panic!("thermal coefficients for {name}: {error}"));
        let (cp_molar, enthalpy, _) = subs
            .calculate_thermo_properties(name, THERMO_TEMPERATURE_K)
            .unwrap_or_else(|error| panic!("thermochemistry for {name}: {error}"));
        let (molar_mass_g_per_mol, _) = calculate_molar_mass(name.clone(), None)
            .unwrap_or_else(|error| panic!("molar mass for {name}: {error}"));
        let source = subs
            .get_search_result(name, WhatIsFound::Thermo)
            .expect("resolved thermo search result")
            .library()
            .to_string();
        thermo.insert(name.clone(), (cp_molar, enthalpy));
        molar_masses.insert(name.clone(), molar_mass_g_per_mol / 1000.0);
        provenance.insert(name.clone(), source);
    }
    (thermo, molar_masses, provenance)
}

fn reaction_heat_release_j_per_mol(equation: &str, thermo: &HashMap<String, (f64, f64)>) -> f64 {
    let (reactants, products) = parse_equation(equation);
    let reactant_enthalpy: f64 = reactants
        .iter()
        .map(|(name, coefficient)| coefficient * thermo[name].1)
        .sum();
    let product_enthalpy: f64 = products
        .iter()
        .map(|(name, coefficient)| coefficient * thermo[name].1)
        .sum();
    reactant_enthalpy - product_enthalpy
}

pub fn build_reactor() -> SimpleReactorTask {
    let (reactions, species) = load_reactions();
    let (thermo, molar_masses, _) = resolved_thermochemistry(&species);
    let reaction_rows: Vec<FastElemReact> = reactions
        .iter()
        .map(|reaction| FastElemReact {
            eq: reaction.equation.clone(),
            A: reaction.a_si,
            n: reaction.temperature_exponent,
            E: reaction.activation_energy_j_per_mol,
            Q: reaction_heat_release_j_per_mol(&reaction.equation, &thermo),
        })
        .collect();

    let mut mole_fractions: HashMap<String, f64> =
        species.iter().cloned().map(|name| (name, 0.0)).collect();
    for (name, fraction) in INLET_MOLE_FRACTIONS {
        *mole_fractions
            .get_mut(name)
            .unwrap_or_else(|| panic!("inlet species {name} is not in the mechanism")) = fraction;
    }
    let mole_fraction_sum: f64 = mole_fractions.values().sum();
    assert!((mole_fraction_sum - 1.0).abs() < 1e-12);
    let mean_molar_mass_kg_per_mol: f64 = species
        .iter()
        .map(|name| mole_fractions[name] * molar_masses[name])
        .sum();
    let inlet_mass_fractions: HashMap<String, f64> = species
        .iter()
        .map(|name| {
            (
                name.clone(),
                mole_fractions[name] * molar_masses[name] / mean_molar_mass_kg_per_mol,
            )
        })
        .collect();
    let cp_mass_j_per_kg_k: f64 = species
        .iter()
        .map(|name| inlet_mass_fractions[name] * thermo[name].0 / molar_masses[name])
        .sum();
    assert!(cp_mass_j_per_kg_k.is_finite() && cp_mass_j_per_kg_k > 0.0);
    let pressure_pa =
        GAS_DENSITY_KG_PER_M3 * R_G * THERMO_TEMPERATURE_K / mean_molar_mass_kg_per_mol;

    let mut boundary = HashMap::from([
        ("T".to_string(), THERMO_TEMPERATURE_K),
        ("q".to_string(), 0.0),
    ]);
    let mut diffusion = HashMap::new();
    for name in &species {
        boundary.insert(name.clone(), inlet_mass_fractions[name]);
        boundary.insert(format!("{name}_flux"), 0.0);
        diffusion.insert(name.clone(), SPECIES_DIFFUSIVITY_M2_PER_S);
    }

    let mut reactor = SimpleReactorTask::new();
    reactor
        .fast_react_set(reaction_rows)
        .expect("Cantera elementary subset should parse through the BVP reaction facade");
    reactor
        .kindata
        .set_species_universe(species.clone())
        .expect("explicit species universe should retain the inert diluent");
    reactor.set_parameters(
        reactor.thermal_effects.clone(),
        pressure_pa,
        THERMO_TEMPERATURE_K,
        cp_mass_j_per_kg_k,
        boundary,
        THERMAL_CONDUCTIVITY_W_PER_M_K,
        diffusion,
        MASS_FLUX_KG_PER_M2_S,
        ScalingConfig::new(THERMO_TEMPERATURE_K, REACTOR_LENGTH_M, TEMPERATURE_SCALE_K),
    );
    reactor.M = mean_molar_mass_kg_per_mol;
    reactor
}

pub fn inlet_initial_guess(reactor: &SimpleReactorTask, nodes: usize) -> DMatrix<f64> {
    let substances = reactor.kindata.substances();
    let state_count = reactor.solver.unknowns.len();
    DMatrix::from_fn(state_count, nodes, |row, _| {
        if row < 2 || row % 2 == 1 {
            0.0
        } else {
            let species_index = (row - 2) / 2;
            let name = &substances[species_index];
            reactor.boundary_condition[name]
        }
    })
}
