use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::Kinetics::kinetics_lib_api::{cached_reactbase, cached_reaction_db};
use crate::Kinetics::stoichiometry_analyzer::clean_off_DUP;
use serde_json::Value;
use std::collections::{HashMap, HashSet};

fn vec_string_sorted(vec_strings: &[String]) -> KineticsResult<Vec<String>> {
    let mut vec_ints: Vec<i32> = vec_strings
        .iter()
        .map(|s| {
            s.parse::<i32>().map_err(|_| {
                KineticsError::InvalidReactionData(format!("invalid mechanism reaction id `{}`", s))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    vec_ints.sort();
    Ok(vec_ints.iter().map(|i| i.to_string()).collect())
}

/// Walk the reaction graph using already loaded in-memory JSON data.
///
/// The input and working sets stay owned, so we do not clone the full database
/// or rebuild a borrowed mirror of it on every call.
pub fn parse_database(
    db_object: &HashMap<String, HashMap<String, HashSet<String>>>,
    mut mechanism: HashSet<String>,
    mut reactants: HashSet<String>,
) -> (HashSet<String>, HashSet<String>) {
    let mut found_reactants: HashSet<String> = HashSet::new();
    let mut found_reactions: HashSet<String> = HashSet::new();
    let keys_to_replace = HashSet::from([
        String::from("M"),
        String::from("M)"),
        String::from("(M"),
        String::from("(M)"),
    ]);

    loop {
        for (r_id, reaction_data) in db_object
            .iter()
            .filter(|(r_id, _)| !mechanism.contains(r_id.as_str()))
        {
            let Some(reactants_db) = reaction_data.get("reagents") else {
                continue;
            };

            let reagents_ok = reactants_db
                .iter()
                .filter(|subs| !keys_to_replace.contains(*subs))
                .all(|subs| reactants.contains(subs));

            if reagents_ok {
                if let Some(products_db) = reaction_data.get("products") {
                    found_reactants.extend(products_db.iter().cloned());
                }
                found_reactions.insert(r_id.clone());
            }
        }

        if found_reactions.is_empty() {
            return (reactants, mechanism);
        }

        mechanism.extend(found_reactions.drain());
        reactants.extend(found_reactants.drain());
    }
}

/// Search the selected mechanism and return discovered reaction ids,
/// cleaned substances, reaction equations, and reaction data payloads.
pub fn mechfinder(
    big_mech: &str,
    vec: Vec<&str>,
) -> KineticsResult<(Vec<String>, Vec<String>, Vec<String>, Vec<Value>)> {
    let reactants: HashSet<String> = vec.iter().map(|s| (*s).to_string()).collect();
    let mechanism: HashSet<String> = HashSet::new();

    let reactlibrary = cached_reactbase()?;
    let reactlibrary = reactlibrary
        .get(big_mech)
        .ok_or_else(|| KineticsError::MissingLibrary(big_mech.to_string()))?;

    let reaction_db = cached_reaction_db()?;
    let reaction_db = reaction_db
        .get(big_mech)
        .ok_or_else(|| KineticsError::MissingLibrary(big_mech.to_string()))?;
    let (mut reactants, mechanism) = parse_database(reaction_db, mechanism, reactants);

    let keys_to_replace = HashSet::from(["M", "M)", "(M", "(M)", "V", "(V", "V)"]);
    reactants.retain(|subs| !keys_to_replace.contains(subs.as_str()));

    let mut reactants: Vec<String> = reactants
        .into_iter()
        .map(|mut item| {
            clean_off_DUP(&mut item);
            item
        })
        .collect();
    reactants.sort();

    let mechanism: Vec<String> = mechanism.into_iter().collect();
    let mechanism = vec_string_sorted(&mechanism)?;

    let mut vec_of_reactions: Vec<String> = Vec::with_capacity(mechanism.len());
    let mut vec_of_reaction_value: Vec<Value> = Vec::with_capacity(mechanism.len());
    for react_num in &mechanism {
        let data_for_react = reactlibrary
            .get(react_num)
            .ok_or_else(|| KineticsError::MissingReaction(react_num.clone()))?;
        vec_of_reaction_value.push(data_for_react.to_owned());
        let json_string = serde_json::to_string(data_for_react.get("eq").ok_or_else(|| {
            KineticsError::InvalidReactionData(format!(
                "reaction `{}` is missing field `eq`",
                react_num
            ))
        })?)?;
        vec_of_reactions.push(json_string);
    }
    Ok((
        mechanism,
        reactants,
        vec_of_reactions,
        vec_of_reaction_value,
    ))
}
