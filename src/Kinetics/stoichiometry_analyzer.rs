//! # Stoichiometry Analyzer Module
//!
//! ## Purpose
//! This module analyzes chemical reaction equations and generates comprehensive stoichiometric data.
//! It parses reaction strings, extracts substances, and creates matrices for kinetic modeling.
//! Handles database parsing artifacts and supports empirical reactions with custom power coefficients.
//!
//! ## Main Data Structures
//! - `StoichAnalyzer`: Core structure containing all stoichiometric analysis results
//!   - `reactions`: Vector of reaction equation strings
//!   - `substances`: Vector of unique substance names found in reactions
//!   - `stecheo_matrx`: Stoichiometric coefficient matrix (reactions × substances)
//!   - `stecheo_reags`/`stecheo_prods`: Separate matrices for reactants and products
//!   - `G_matrix_reag`/`G_matrix_prod`: Power coefficient matrices for kinetic equations
//!   - `matrix_of_elements`: Elemental composition matrix
//!   - `vec_of_molmasses`: Molecular masses of substances
//!
//! ## Key Logic Implementation
//! 1. **Artifact Cleaning**: Removes duplicate-ID suffixes and third-body markers without erasing formula parentheses or phase identity
//! 2. **Reaction Parsing**: Handles supported arrows, charged species, and balanced parenthesized names
//! 3. **Substance Extraction**: Parses positive decimal/scientific coefficients and explicit kinetic powers (A**0.3)
//! 4. **Matrix Construction**: Builds multiple matrices for different aspects of kinetic modeling
//! 5. **Third Body Handling**: Recognizes and filters out collision partners marked as 'M'
//!
//! ## Usage Pattern
//! ```rust
//! use KiThe::Kinetics::stoichiometry_analyzer::StoichAnalyzer;
//! let mut analyzer = StoichAnalyzer::new();
//! analyzer.reactions = vec!["A=2B".to_string(), "B->A + 3C".to_string()];
//! analyzer.search_substances().unwrap();  // Find all unique substances
//! analyzer.analyse_reactions().unwrap();  // Generate stoichiometric matrices
//! analyzer.create_matrix_of_elements().unwrap();  // Generate elemental composition
//! ```
//!
//! ## Interesting Features
//! - **Empirical Reaction Support**: Handles power coefficients different from stoichiometric (A**0.5)
//! - **Database Artifact Cleaning**: Removes '_DUP', '_dup' suffixes and collision partner markers
//! - **Flexible Arrow Recognition**: Supports =, ->, =>, <=>, and variations with spaces
//! - **Third Body Recognition**: Automatically filters out 'M' and '(M)' collision partners
//! - **Dual Matrix System**: Separate tracking of stoichiometric vs kinetic power coefficients
//! - **Integrated Molecular Data**: Combines with molmass module for complete substance characterization

use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::Kinetics::molmass::{
    calculate_molar_mass_of_vector_of_subs, create_elem_composition_matrix,
};
use nalgebra::DMatrix;
use regex::Regex;
/// ru
/// Модуль берет на вход  вектор уравнений реации, заданных в виде String и выдает следующие данные:
/// 1) стехиометрическую матрицу заданную в виде вектора векторов
/// 2) вектор веществ
/// 3) вектор векторов стехиометрических коэффициентов реагентов в каждой реакции
/// 4) то же для продуктов
/// В процессе производится избавление от артефактов парсинга уравнений из баз данных
/// Примечание:
/// 1) уравнения реакции содержат последние символы '_dup' или '_DUP' (то есть дублирующие) - эти
/// символы удаляются как артефакты парсинга
/// 2) для работы с эмпирическими реакциями у которых
///  код возвращает следующие структуры данных: матрица стехеометрических коэффициентов,
///  матрица коэффициентов прямых реакций и матрица коэффициентов обратных реакций, матрица степеней концентраций для кинетической функции, G_matrix,
/// как правило степени концентраций в кинетической функции совпадают со стехеометрическими коэффициентами веществ в реакциии, однако,
/// для эмпирических реакций они могут и отличаться от стехеометрических коэффициентов.
///  Предусмотрена вохможность прямого указания этих коэффициентов в уравнениях реакции в виде "степени"
/// после формулы вещества (например A**0.3 ) такие степени записываются в  G_matrix вместо стехеометрических коэффициентов
/// в противном случае записывается стехеометрический коэффициент
/// ----------------------------------------------------------------
/// eng
/// The module takes as input a vector of reaction equations specified as a vector of String and produces the following data:
/// 1) a stoichiometric matrix specified as a vector of vectors
/// 2) a vector of substances
/// 3) a vector of vectors of stoichiometric coefficients of reactants in each reaction
/// 4) the same for products
/// The process involves getting rid of artifacts of parsing equations from databases
///
/// Note:
/// 1) reaction equations contain the last characters '_dup' or '_DUP' (that is, duplicate ones) - these
/// characters are removed as parsing artifacts
/// 2) for working with empirical reactions in which  the code returns the following data structures: matrix of stoicheometric coefficients,
/// matrix of coefficients of direct reactions and matrix of coefficients of reverse reactions, matrix of degrees of concentration for the
/// kinetic function, G_matrix. As a rule, the degrees of concentration in the kinetic function coincide with the stoicheometric coefficients of
/// the substances in the reaction; however, for empirical reactions they may differ from the stoicheometric coefficients.
/// It is possible to directly indicate these coefficients in the reaction equations in the form of a “degree” after the formula of a substance
/// (for example A**0.3), such degrees are written in G_matrix instead of stoicheometric coefficients otherwise, the stoicheometric coefficient is written
/// # Examples
/// reaction data with  artifacts of parsing equations from databases
/// ```
/// use KiThe::Kinetics::stoichiometry_analyzer::StoichAnalyzer  ;
/// let mut ReactionAnalyzer_instance = StoichAnalyzer  ::new();
/// let reactions_: Vec<&str> = vec!["A=2BM)", "B->A + 3C_DUP", "2B+A=D"];
/// let reaction = reactions_.iter().map(|s| s.to_string()).collect();
/// ReactionAnalyzer_instance.reactions = reaction;
/// ReactionAnalyzer_instance.search_substances().unwrap();
/// println!("substances: {:?}", ReactionAnalyzer_instance.substances);
/// let substancses_ = vec!["A", "B", "C", "D"];
/// let substancses = substancses_.iter().map(|s| s.to_string()).collect();
/// ReactionAnalyzer_instance.substances = substancses ;
/// ReactionAnalyzer_instance.analyse_reactions().unwrap();
/// println!("{:?}", ReactionAnalyzer_instance);
/// ```
use std::collections::{HashMap, HashSet};

fn reaction_syntax(input: &str, reason: impl Into<String>) -> KineticsError {
    KineticsError::ReactionSyntax {
        input: input.to_string(),
        reason: reason.into(),
    }
}

fn split_reaction(reaction: &str) -> KineticsResult<(&str, &str)> {
    let arrow_re =
        Regex::new(r"<=\s*>|<\s*=>|<=>|=>|->|=").map_err(|source| KineticsError::RegexCompile {
            context: "reaction arrow",
            source,
        })?;
    let mut arrows = arrow_re.find_iter(reaction);
    let arrow = arrows
        .next()
        .ok_or_else(|| reaction_syntax(reaction, "no supported reaction arrow"))?;
    if arrows.next().is_some() {
        return Err(reaction_syntax(reaction, "more than one reaction arrow"));
    }
    let left = reaction[..arrow.start()].trim();
    let right = reaction[arrow.end()..].trim();
    if left.is_empty() || right.is_empty() {
        return Err(reaction_syntax(
            reaction,
            "reaction sides must not be empty",
        ));
    }
    Ok((left, right))
}

struct ParsedParticipant {
    name: String,
    coefficient: f64,
    kinetic_order: f64,
}

struct ParsedReaction {
    equation: String,
    reactants: Vec<ParsedParticipant>,
    products: Vec<ParsedParticipant>,
}

fn parse_reaction(input: &str) -> KineticsResult<ParsedReaction> {
    let equation = input.replace("**", "^").trim().to_string();
    let mut cleaned = equation.clone();
    let equation_for_parse = clean_off_DUP(&mut cleaned);
    let (left, right) = split_reaction(equation_for_parse)?;
    let parse_side = |side: &str, side_name: &str| -> KineticsResult<Vec<ParsedParticipant>> {
        let (names, coefficients, orders) = analyse_substances(side)
            .map_err(|error| error.with_context(format!("{side_name} side of `{equation}`")))?;
        Ok(names
            .into_iter()
            .zip(coefficients)
            .zip(orders)
            .map(|((name, coefficient), kinetic_order)| ParsedParticipant {
                name,
                coefficient,
                kinetic_order,
            })
            .collect())
    };

    let reactants = parse_side(left, "reactant")?;
    let products = parse_side(right, "product")?;
    Ok(ParsedReaction {
        equation,
        reactants,
        products,
    })
}

fn parse_reaction_rows(
    reactions: &[String],
    reaction_contexts: Option<&[String]>,
) -> KineticsResult<Vec<ParsedReaction>> {
    if let Some(contexts) = reaction_contexts {
        if contexts.len() != reactions.len() {
            return Err(KineticsError::LengthMismatch {
                context: "reaction equations and diagnostic identities",
                left: reactions.len(),
                right: contexts.len(),
            });
        }
    }
    reactions
        .iter()
        .enumerate()
        .map(|(index, equation)| {
            parse_reaction(equation).map_err(|error| {
                let location = reaction_contexts
                    .and_then(|contexts| contexts.get(index))
                    .map(|context| format!(" ({context})"))
                    .unwrap_or_default();
                error.with_context(format!("reaction row {index}{location}: `{equation}`"))
            })
        })
        .collect()
}

/*
        #  строковая переменная на которые разделяет по знаку "+" половину уравнения относящююся к
        # реагентам/продуктам может содержать '(', ')', а также 'M' - символ третьей частицы которая не относится к
        # продуктам и реагентам
*/
fn clean_third_body_marker(species: &str) -> Option<&str> {
    if matches!(species, "M" | "m" | "(M)" | "(+M)" | "(+m)") {
        return None;
    }
    let species = species
        .strip_suffix("(+M)")
        .or_else(|| species.strip_suffix("(+m)"))
        .or_else(|| species.strip_suffix("(M)"))
        .or_else(|| species.strip_suffix("M)"))
        .unwrap_or(species)
        .trim();
    (!species.is_empty()).then_some(species)
}

fn validate_species_parentheses(species: &str, source: &str) -> KineticsResult<()> {
    let mut depth = 0usize;
    for character in species.chars() {
        match character {
            '(' => depth += 1,
            ')' if depth == 0 => {
                return Err(reaction_syntax(source, "unmatched closing parenthesis"));
            }
            ')' => depth -= 1,
            _ => {}
        }
    }
    if depth != 0 {
        return Err(reaction_syntax(
            source,
            "unclosed parenthesis in species name",
        ));
    }
    Ok(())
}

fn has_explicit_charge_magnitude(species: &str) -> bool {
    let Some((_, suffix)) = species.rsplit_once('^') else {
        return false;
    };
    let magnitude = suffix
        .strip_suffix('+')
        .or_else(|| suffix.strip_suffix('-'));
    magnitude.is_some_and(|value| !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit()))
}

fn split_participants(side: &str) -> KineticsResult<Vec<&str>> {
    let characters: Vec<(usize, char)> = side.char_indices().collect();
    let mut pieces = Vec::new();
    let mut start = 0;
    let mut depth = 0usize;
    for (index, (byte, character)) in characters.iter().copied().enumerate() {
        match character {
            '(' => depth += 1,
            ')' if depth == 0
                && index + 1 == characters.len()
                && characters[index.saturating_sub(1)].1 == 'M' => {}
            ')' => {
                if depth == 0 {
                    return Err(reaction_syntax(side, "unmatched closing parenthesis"));
                }
                depth -= 1;
            }
            '+' if depth == 0 => {
                let next = characters[index + 1..]
                    .iter()
                    .find(|(_, ch)| !ch.is_whitespace())
                    .map(|(_, ch)| *ch);
                let previous_is_species = side[..byte]
                    .chars()
                    .next_back()
                    .is_some_and(|ch| !ch.is_whitespace());
                let is_charge = matches!(next, None | Some('=') | Some('>'))
                    || (next == Some('+') && previous_is_species);
                if !is_charge {
                    pieces.push(&side[start..byte]);
                    start = byte + character.len_utf8();
                }
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(reaction_syntax(side, "unclosed parenthesis"));
    }
    pieces.push(&side[start..]);
    Ok(pieces)
}

pub fn analyse_substances(
    half_reaction: &str,
) -> KineticsResult<(Vec<String>, Vec<f64>, Vec<f64>)> {
    // Participant separators are plus signs outside parentheses; a terminal
    // plus on an ion is retained as part of the species name.
    let mut s_list = Vec::new();
    let mut g_list = Vec::new();
    let mut subs: Vec<String> = Vec::new();

    let re_coeff =
        Regex::new(r"^((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\*?\s*").map_err(|source| {
            KineticsError::RegexCompile {
                context: "stoichiometric coefficient",
                source,
            }
        })?;
    //let re_power = Regex::new(r"\^\ *(\d+(\.\d*)?)").unwrap();
    let re_power =
        Regex::new(r"\^((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$").map_err(|source| {
            KineticsError::RegexCompile {
                context: "kinetic order",
                source,
            }
        })?;

    for s in split_participants(half_reaction)? {
        let mut s = s.trim();
        if s.is_empty() {
            return Err(reaction_syntax(half_reaction, "empty participant"));
        }
        let mut stec_coeff: f64 = 1.0;
        let _end_of_stoichiometric: usize = 0;
        //  captures.get(1) is used to access the first capture group in the regular expression match.
        if let Some(captures) = re_coeff.captures(s) {
            // refers to the first capture group in the regular expression match, example: Match { start: 0, end: 1, string: "2" }
            let regmatch = captures
                .get(0)
                .ok_or_else(|| reaction_syntax(half_reaction, "invalid stoichiometric prefix"))?;
            let start_of_stoichiometric = regmatch.start();
            let end_of_stoichiometric = regmatch.end();
            // стехиометрический коэффициент должен быть перед формулой вещества, т.е. номер его позиции должен начинаться с нуля
            if start_of_stoichiometric == 0 {
                s = &s[end_of_stoichiometric..];
                stec_coeff = captures
                    .get(1)
                    .ok_or_else(|| {
                        reaction_syntax(half_reaction, "missing stoichiometric coefficient")
                    })?
                    .as_str()
                    .parse()
                    .map_err(|_| {
                        reaction_syntax(half_reaction, "invalid stoichiometric coefficient")
                    })?;
                if !stec_coeff.is_finite() || stec_coeff <= 0.0 {
                    return Err(reaction_syntax(
                        half_reaction,
                        "stoichiometric coefficient must be finite and positive",
                    ));
                }
            }
        }
        // Default the kinetic power to the stoichiometric coefficient until
        // we see an explicit power override in the reaction text.
        let mut power_coeff: f64 = stec_coeff;

        let charged_with_magnitude = has_explicit_charge_magnitude(s);
        if let Some(captures) = re_power.captures(s).filter(|_| !charged_with_magnitude) {
            let regmatch = captures
                .get(0)
                .ok_or_else(|| reaction_syntax(half_reaction, "invalid kinetic order"))?;
            power_coeff = captures
                .get(1)
                .ok_or_else(|| reaction_syntax(half_reaction, "missing kinetic order"))?
                .as_str()
                .parse()
                .map_err(|_| reaction_syntax(half_reaction, "invalid kinetic order"))?;
            if !power_coeff.is_finite() || power_coeff <= 0.0 {
                return Err(reaction_syntax(
                    half_reaction,
                    "kinetic order must be finite and positive",
                ));
            }
            let start = regmatch.start();
            #[allow(unused_variables)]
            let end = regmatch.end();
            s = &s[..start];
        }

        s = s.trim();
        if s.starts_with('/') {
            return Err(reaction_syntax(
                half_reaction,
                "rational stoichiometric coefficients are unsupported; use decimal notation",
            ));
        }
        if s.contains('^') && !charged_with_magnitude {
            return Err(reaction_syntax(half_reaction, "invalid kinetic order"));
        }
        let Some(s) = clean_third_body_marker(s) else {
            continue;
        };
        validate_species_parentheses(s, half_reaction)?;
        let s = s.to_string();
        s_list.push(stec_coeff);
        g_list.push(power_coeff);
        subs.push(s);
    }

    Ok((subs, s_list, g_list))
}

// #  некоторые записи уравнения реакции содержат последние символы '_dup' или '_DUP' (то есть дублирующие)
// избавляемся от них
pub fn clean_off_DUP(item: &mut String) -> &String {
    if item.ends_with("_dup") || item.ends_with("_DUP") {
        *item = item.replace("_dup", "").replace("_DUP", "");
        return item;
    } else {
        item
    }
}

/// This struct is used to represent a reaction.
#[derive(Debug, Clone)]
pub struct StoichAnalyzer {
    pub reactions: Vec<String>, // вектор реакций/ a vector of reactions
    pub groups: Option<HashMap<String, HashMap<String, usize>>>, // Chemical formulae may contain spectial names for chemical groupls i.e. groups of atoms, e.g. Me (methyl) group, which is converted into {"C":1, "H":3}
    pub substances: Vec<String>, // вектор веществ/ a vector of substances
    pub stecheo_matrx: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов в каждой реакции/ a vector of vectors of stoichiometric coefficients in each reaction
    pub stecheo_reags: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов реагентов в каждой реакции/ a vector of vectors of stoichiometric coefficients of reactants in each reaction
    pub stecheo_prods: Vec<Vec<f64>>, // вектор векторов стехиометрических коэффициентов продуктов в каждой реакции/ a vector of vectors of stoichiometric coefficients of products in each reaction
    pub G_matrix_reag: Vec<Vec<f64>>, // matrix of powers of concentrations for constructng kinetic equation for forward reaction
    pub G_matrix_prod: Vec<Vec<f64>>, // matrix of powers of concentrations for constructng kinetic equation for reverse reaction
    pub matrix_of_elements: Option<DMatrix<f64>>, // matrix of elemental composition
    pub vec_of_molmasses: Option<Vec<f64>>, // vector of molecular masses
    pub unique_vec_of_elems: Option<Vec<String>>,
    composition_source: Option<(Vec<String>, Option<HashMap<String, HashMap<String, usize>>>)>,
}

impl StoichAnalyzer {
    pub fn new() -> Self {
        StoichAnalyzer {
            reactions: Vec::new(),
            groups: None,
            substances: Vec::new(),
            stecheo_matrx: Vec::new(),
            stecheo_reags: Vec::new(),
            stecheo_prods: Vec::new(),
            G_matrix_reag: Vec::new(),
            G_matrix_prod: Vec::new(),
            matrix_of_elements: None,
            vec_of_molmasses: None,
            unique_vec_of_elems: None,
            composition_source: None,
        }
    }

    pub fn analyse_reactions(&mut self) -> KineticsResult<()> {
        self.analyse_reactions_with_context(None)
    }

    pub(crate) fn analyse_reactions_with_context(
        &mut self,
        reaction_contexts: Option<&[String]>,
    ) -> KineticsResult<()> {
        self.invalidate_composition_if_source_changed();
        let parsed_reactions = parse_reaction_rows(&self.reactions, reaction_contexts)?;
        let rows = parsed_reactions.len();
        let columns = self.substances.len();
        let species_index: HashMap<&str, usize> = self
            .substances
            .iter()
            .enumerate()
            .map(|(index, name)| (name.as_str(), index))
            .collect();
        let mut stoich = vec![vec![0.0; columns]; rows];
        let mut reactants = stoich.clone();
        let mut products = stoich.clone();
        let mut reactant_orders = stoich.clone();
        let mut product_orders = stoich.clone();
        for (row, reaction) in parsed_reactions.iter().enumerate() {
            for (participants, is_reactant) in
                [(&reaction.reactants, true), (&reaction.products, false)]
            {
                for participant in participants {
                    let column = species_index
                        .get(participant.name.as_str())
                        .ok_or_else(|| KineticsError::UnknownReactionSpecies {
                            reaction: reaction.equation.clone(),
                            species: participant.name.clone(),
                        })?;
                    if is_reactant {
                        stoich[row][*column] -= participant.coefficient;
                        reactants[row][*column] += participant.coefficient;
                        reactant_orders[row][*column] += participant.kinetic_order;
                    } else {
                        stoich[row][*column] += participant.coefficient;
                        products[row][*column] += participant.coefficient;
                        product_orders[row][*column] += participant.kinetic_order;
                    }
                }
            }
        } // end of for reaction in &self.reactions {
        self.reactions = parsed_reactions
            .into_iter()
            .map(|reaction| reaction.equation)
            .collect();
        self.stecheo_matrx = stoich;
        self.stecheo_reags = reactants;
        self.stecheo_prods = products;
        self.G_matrix_reag = reactant_orders;
        self.G_matrix_prod = product_orders;
        Ok(())
    }

    pub fn search_substances(&mut self) -> KineticsResult<()> {
        self.search_substances_with_context(None)
    }

    pub(crate) fn search_substances_with_context(
        &mut self,
        reaction_contexts: Option<&[String]>,
    ) -> KineticsResult<()> {
        let parsed_reactions = parse_reaction_rows(&self.reactions, reaction_contexts)?;
        let mut found_substances: Vec<String> = Vec::new();
        let mut seen = HashSet::new();
        for reaction in &parsed_reactions {
            for participant in reaction.reactants.iter().chain(&reaction.products) {
                if seen.insert(participant.name.clone()) {
                    found_substances.push(participant.name.clone());
                }
            }
        } // end of for reaction in &self.reactions {
        if self.substances != found_substances {
            self.invalidate_composition_data();
        }
        self.reactions = parsed_reactions
            .into_iter()
            .map(|reaction| reaction.equation)
            .collect();
        self.substances = found_substances;
        Ok(())
    }
    pub fn create_matrix_of_elements(&mut self) -> KineticsResult<()> {
        self.invalidate_composition_if_source_changed();
        let formulae: Vec<String> = self
            .substances
            .iter()
            .map(|species| composition_formula(species))
            .collect();
        let formula_refs: Vec<&str> = formulae.iter().map(String::as_str).collect();
        let computed_masses =
            calculate_molar_mass_of_vector_of_subs(formula_refs.clone(), self.groups.clone())?;
        let (matrix, vec_of_elems) =
            create_elem_composition_matrix(formula_refs, self.groups.clone())?;
        if let Some(masses) = &self.vec_of_molmasses {
            if masses.len() != self.substances.len() {
                return Err(KineticsError::LengthMismatch {
                    context: "species and explicitly supplied molar masses",
                    left: self.substances.len(),
                    right: masses.len(),
                });
            }
        }
        let masses = self.vec_of_molmasses.clone().unwrap_or(computed_masses);
        self.matrix_of_elements = Some(matrix);
        self.vec_of_molmasses = Some(masses);
        self.unique_vec_of_elems = Some(vec_of_elems);
        self.composition_source = Some((self.substances.clone(), self.groups.clone()));
        Ok(())
    }

    /// Replace the species universe and invalidate all composition-derived data.
    pub fn set_substances(&mut self, substances: Vec<String>) {
        if self.substances != substances {
            self.substances = substances;
            self.invalidate_composition_data();
        }
    }

    fn invalidate_composition_if_source_changed(&mut self) {
        if self
            .composition_source
            .as_ref()
            .is_some_and(|(species, groups)| species != &self.substances || groups != &self.groups)
        {
            self.invalidate_composition_data();
        }
    }

    fn invalidate_composition_data(&mut self) {
        self.matrix_of_elements = None;
        self.vec_of_molmasses = None;
        self.unique_vec_of_elems = None;
        self.composition_source = None;
    }
} // end of impl ReactionAnalyzer

fn composition_formula(species: &str) -> String {
    let phase_stripped = [
        "(aq)", "(g)", "(G)", "(l)", "(L)", "(s)", "(S)", "(c)", "(C)", "(cr)",
    ]
    .iter()
    .find_map(|suffix| species.strip_suffix(suffix))
    .unwrap_or(species);
    if has_explicit_charge_magnitude(phase_stripped) {
        return phase_stripped
            .rsplit_once('^')
            .map(|(formula, _)| formula)
            .unwrap_or(phase_stripped)
            .to_string();
    }
    phase_stripped
        .strip_suffix('+')
        .or_else(|| phase_stripped.strip_suffix('-'))
        .unwrap_or(phase_stripped)
        .to_string()
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TESTS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analyse_substances() {
        let half_reaction = "5H2O + 10O2";
        let (subs, s_list, g_list) = analyse_substances(half_reaction).unwrap();
        assert_eq!(subs, vec!["H2O".to_string(), "O2".to_string()]);
        assert_eq!(s_list, vec![5.0, 10.0]);
        assert_eq!(g_list, vec![5.0, 10.0]);
    }

    #[test]
    fn test_analyse_substances_prefers_explicit_powers() {
        let half_reaction = "2H2O^3 + 4O2";
        let (subs, s_list, g_list) = analyse_substances(half_reaction).unwrap();
        assert_eq!(subs, vec!["H2O".to_string(), "O2".to_string()]);
        assert_eq!(s_list, vec![2.0, 4.0]);
        assert_eq!(g_list, vec![3.0, 4.0]);
    }

    #[test]
    fn test_analyse_substances_with_exotic_names() {
        let mut ReactionAnalyzer_instance = StoichAnalyzer::new();
        let reactions_: Vec<&str> = vec!["HMX=>7HMXprod"];
        let reaction = reactions_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.reactions = reaction;
        // ReactionAnalyzer_instance.search_substances();
        // assert_eq!(ReactionAnalyzer_instance.substances, vec!["A".to_string(), "B".to_string(), "C".to_string()]);
        let substancses_ = vec!["HMX", "HMXprod"];
        let substancses = substancses_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.substances = substancses;
        ReactionAnalyzer_instance.analyse_reactions().unwrap();
        let stecheo_matrx = ReactionAnalyzer_instance.stecheo_matrx;
        let result = [[-1.0, 7.0]];
        let result: Vec<Vec<f64>> = result.iter().map(|row| row.to_vec()).collect();
        assert_eq!(stecheo_matrx, result);
    }
    #[test]
    fn test_clean_off_DUP() {
        let mut item = "B->A + 3C_DUP".to_string();
        let cleaned_item = clean_off_DUP(&mut item);
        assert_eq!(cleaned_item, "B->A + 3C");
    }

    #[test]
    fn test_search_substances() {
        let mut ReactionAnalyzer_instance = StoichAnalyzer::new();
        let reactions_: Vec<&str> = vec!["A=2BM)", "B->A + 3C_DUP"];
        let reaction = reactions_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.reactions = reaction;
        ReactionAnalyzer_instance.search_substances().unwrap();
        assert_eq!(
            ReactionAnalyzer_instance.substances,
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
    }
    #[test]
    fn test_reaction_analyzer() {
        let mut ReactionAnalyzer_instance = StoichAnalyzer::new();
        let reactions_: Vec<&str> = vec!["A=2BM)", "B=>A + 3C_DUP", "2B+A=D**0.5"];
        let reaction = reactions_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.reactions = reaction;
        // ReactionAnalyzer_instance.search_substances();
        // assert_eq!(ReactionAnalyzer_instance.substances, vec!["A".to_string(), "B".to_string(), "C".to_string()]);
        let substancses_ = vec!["A", "B", "C", "D"];
        let substancses = substancses_.iter().map(|s| s.to_string()).collect();
        ReactionAnalyzer_instance.substances = substancses;
        ReactionAnalyzer_instance.analyse_reactions().unwrap();
        let stecheo_matrx = ReactionAnalyzer_instance.stecheo_matrx;
        let result = [
            [-1.0, 2.0, 0.0, 0.0],
            [1.0, -1.0, 3.0, 0.0],
            [-1.0, -2.0, 0.0, 1.0],
        ];
        let result: Vec<Vec<f64>> = result.iter().map(|row| row.to_vec()).collect();
        assert_eq!(stecheo_matrx, result);
    }

    #[test]
    fn repeated_analysis_and_duplicate_reaction_rows_are_stable() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = vec!["A=>B".into(), "A=>B".into()];
        analyzer.substances = vec!["A".into(), "B".into()];

        analyzer.analyse_reactions().unwrap();
        let first = analyzer.stecheo_matrx.clone();
        analyzer.analyse_reactions().unwrap();

        assert_eq!(analyzer.stecheo_matrx, vec![vec![-1.0, 1.0]; 2]);
        assert_eq!(analyzer.stecheo_matrx, first);
        assert_eq!(analyzer.stecheo_reags, vec![vec![1.0, 0.0]; 2]);
        assert_eq!(analyzer.stecheo_prods, vec![vec![0.0, 1.0]; 2]);
    }

    #[test]
    fn repeated_participants_sum_stoichiometry_and_kinetic_orders() {
        let (names, coefficients, orders) =
            analyse_substances("H + H^0.5 + 2H").expect("repeated participants are valid");
        assert_eq!(names, vec!["H", "H", "H"]);
        assert_eq!(coefficients, vec![1.0, 1.0, 2.0]);
        assert_eq!(orders, vec![1.0, 0.5, 2.0]);

        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = vec!["H + H^0.5 + 2H => H2".into()];
        analyzer.substances = vec!["H".into(), "H2".into()];
        analyzer.analyse_reactions().unwrap();
        assert_eq!(analyzer.stecheo_matrx, vec![vec![-4.0, 1.0]]);
        assert_eq!(analyzer.stecheo_reags, vec![vec![4.0, 0.0]]);
        assert_eq!(analyzer.G_matrix_reag, vec![vec![3.5, 0.0]]);
    }

    #[test]
    fn parser_preserves_charges_phases_and_formula_parentheses() {
        let (species, coefficients, _) =
            analyse_substances("H+ + OH-(aq) + Fe^3+ + SO4^2- + 0.5Ca(OH)2(s) + A(+M)")
                .expect("charged, phased, and parenthesized participants are valid");
        assert_eq!(
            species,
            vec!["H+", "OH-(aq)", "Fe^3+", "SO4^2-", "Ca(OH)2(s)", "A"]
        );
        assert_eq!(coefficients, vec![1.0, 1.0, 1.0, 1.0, 0.5, 1.0]);

        let (species, coefficients, _) =
            analyse_substances("1e-3H2 + 2H2O").expect("scientific coefficients are valid");
        assert_eq!(species, vec!["H2", "H2O"]);
        assert_eq!(coefficients, vec![1e-3, 2.0]);
    }

    #[test]
    fn parser_rejects_zero_coefficients_and_unbalanced_species_parentheses() {
        assert!(matches!(
            analyse_substances("0H2"),
            Err(KineticsError::ReactionSyntax { .. })
        ));
        assert!(matches!(
            analyse_substances("1/2H2"),
            Err(KineticsError::ReactionSyntax { .. })
        ));
        assert!(matches!(
            analyse_substances("Ca(OH2"),
            Err(KineticsError::ReactionSyntax { .. })
        ));
    }

    #[test]
    fn reaction_parser_rejects_missing_multiple_and_empty_sides() {
        assert!(matches!(
            split_reaction("A + B"),
            Err(KineticsError::ReactionSyntax { input, reason })
                if input == "A + B" && reason.contains("arrow")
        ));
        assert!(matches!(
            split_reaction("A => B = C"),
            Err(KineticsError::ReactionSyntax { reason, .. })
                if reason.contains("more than one")
        ));
        assert!(matches!(
            split_reaction("=> B"),
            Err(KineticsError::ReactionSyntax { reason, .. })
                if reason.contains("empty")
        ));
        assert!(matches!(
            analyse_substances("A + + B"),
            Err(KineticsError::ReactionSyntax { input, .. }) if input == "A + + B"
        ));
    }

    #[test]
    fn failed_analysis_does_not_publish_partial_matrices_or_drop_unknown_species() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = vec!["A=>B".into(), "B=>C".into()];
        analyzer.substances = vec!["A".into(), "B".into()];
        analyzer.stecheo_matrx = vec![vec![42.0]];
        analyzer.G_matrix_reag = vec![vec![17.0]];

        assert!(matches!(
            analyzer.analyse_reactions(),
            Err(KineticsError::UnknownReactionSpecies { reaction, species })
                if reaction == "B=>C" && species == "C"
        ));
        assert_eq!(analyzer.stecheo_matrx, vec![vec![42.0]]);
        assert_eq!(analyzer.G_matrix_reag, vec![vec![17.0]]);
        assert_eq!(analyzer.reactions, vec!["A=>B", "B=>C"]);
    }

    #[test]
    fn substance_discovery_deduplicates_within_and_across_reactions() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = vec!["A + A => B".into(), "B => A".into()];
        analyzer.search_substances().unwrap();
        assert_eq!(analyzer.substances, vec!["A", "B"]);
    }

    #[test]
    fn matrix_rebuild_tracks_dimension_changes_and_balanced_reference_reaction() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.reactions = vec!["2H2 + O2 => 2H2O".into(), "A => B".into()];
        analyzer.substances = vec![
            "H2".into(),
            "O2".into(),
            "H2O".into(),
            "A".into(),
            "B".into(),
        ];
        analyzer.analyse_reactions().unwrap();
        assert_eq!(analyzer.stecheo_matrx.len(), 2);
        assert!(analyzer.stecheo_matrx.iter().all(|row| row.len() == 5));

        let hydrogen_balance =
            analyzer.stecheo_matrx[0][0] * 2.0 + analyzer.stecheo_matrx[0][2] * 2.0;
        let oxygen_balance = analyzer.stecheo_matrx[0][1] * 2.0 + analyzer.stecheo_matrx[0][2];
        assert_eq!(hydrogen_balance, 0.0);
        assert_eq!(oxygen_balance, 0.0);
        assert_eq!(analyzer.G_matrix_reag[0], analyzer.stecheo_reags[0]);
        assert_eq!(analyzer.G_matrix_prod[0], analyzer.stecheo_prods[0]);

        analyzer.reactions = vec!["A => B".into()];
        analyzer.substances = vec!["A".into(), "B".into()];
        analyzer.analyse_reactions().unwrap();
        assert_eq!(analyzer.stecheo_matrx, vec![vec![-1.0, 1.0]]);
        assert_eq!(analyzer.G_matrix_reag, vec![vec![1.0, 0.0]]);
    }

    #[test]
    fn elemental_and_mass_data_follow_species_and_group_revisions() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.set_substances(vec!["H2".into(), "O2".into()]);
        analyzer.create_matrix_of_elements().unwrap();
        let initial_masses = analyzer.vec_of_molmasses.as_ref().unwrap();
        assert!((initial_masses[0] - 2.016).abs() < 0.001);
        assert!((initial_masses[1] - 31.998).abs() < 0.001);
        assert_eq!(analyzer.matrix_of_elements.as_ref().unwrap().ncols(), 2);

        analyzer.substances = vec!["CO2".into(), "Ca(OH)2".into()];
        analyzer.analyse_reactions().unwrap();
        assert!(analyzer.vec_of_molmasses.is_none());
        assert!(analyzer.matrix_of_elements.is_none());
        analyzer.create_matrix_of_elements().unwrap();
        assert_eq!(analyzer.vec_of_molmasses.as_ref().unwrap().len(), 2);
        assert_eq!(
            analyzer.unique_vec_of_elems.as_ref().unwrap(),
            &["C", "Ca", "H", "O"]
        );

        analyzer.groups = Some(HashMap::from([(
            "Fuel".into(),
            HashMap::from([("C".into(), 1), ("H".into(), 4)]),
        )]));
        analyzer.set_substances(vec!["Fuel".into()]);
        analyzer.create_matrix_of_elements().unwrap();
        assert_eq!(analyzer.unique_vec_of_elems.as_ref().unwrap(), &["C", "H"]);

        analyzer.set_substances(vec!["H2".into()]);
        analyzer.vec_of_molmasses = Some(vec![99.0]);
        analyzer.create_matrix_of_elements().unwrap();
        assert_eq!(analyzer.vec_of_molmasses.as_ref().unwrap(), &[99.0]);
    }

    #[test]
    fn composition_projection_strips_charge_and_phase_but_species_identity_keeps_them() {
        let mut analyzer = StoichAnalyzer::new();
        analyzer.set_substances(vec![
            "H+".into(),
            "OH-(aq)".into(),
            "Fe^3+".into(),
            "SO4^2-".into(),
        ]);
        analyzer.create_matrix_of_elements().unwrap();
        assert_eq!(
            analyzer.substances,
            vec!["H+", "OH-(aq)", "Fe^3+", "SO4^2-"]
        );
        assert_eq!(
            analyzer.unique_vec_of_elems.as_ref().unwrap(),
            &["Fe", "H", "O", "S"]
        );
        assert_eq!(analyzer.vec_of_molmasses.as_ref().unwrap().len(), 4);
    }
}
