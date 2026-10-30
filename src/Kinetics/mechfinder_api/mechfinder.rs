use crate::Kinetics::error::{KineticsError, KineticsResult};
use crate::Kinetics::kinetics_lib_api::cached_catalog_pair;
use crate::Kinetics::stoichiometry_analyzer::clean_off_DUP;
use serde_json::Value;
use std::collections::{HashMap, HashSet};

fn vec_string_sorted(mechanism: &str, ids: &[String]) -> KineticsResult<Vec<String>> {
    let mut sorted = ids.to_vec();
    if let Some(reaction_id) = sorted
        .iter()
        .find(|id| id.is_empty() || !id.bytes().all(|byte| byte.is_ascii_digit()))
    {
        return Err(KineticsError::InvalidMechanismReactionId {
            mechanism: mechanism.to_string(),
            reaction_id: reaction_id.clone(),
        });
    }
    sorted.sort_by(|left, right| {
        let left_significant = left.trim_start_matches('0');
        let right_significant = right.trim_start_matches('0');
        left_significant
            .len()
            .cmp(&right_significant.len())
            .then_with(|| left_significant.cmp(right_significant))
            .then_with(|| left.cmp(right))
    });
    Ok(sorted)
}

fn reaction_equation<'a>(reaction_id: &str, payload: &'a Value) -> KineticsResult<&'a str> {
    payload
        .get("eq")
        .and_then(Value::as_str)
        .ok_or_else(|| KineticsError::InvalidReactionPayload {
            reaction: reaction_id.to_string(),
            field: "eq".to_string(),
            reason: "expected a JSON string".to_string(),
        })
}

/// Walk the directed reaction graph using already loaded in-memory JSON data.
///
/// A reaction is reachable when all non-inert reagents are known; its products
/// are then added to the known species. Empty initial species select no rows.
/// The input and working sets stay owned, so the full database is not cloned.
pub fn parse_database(
    db_object: &HashMap<String, HashMap<String, HashSet<String>>>,
    mut mechanism: HashSet<String>,
    mut reactants: HashSet<String>,
) -> KineticsResult<(HashSet<String>, HashSet<String>)> {
    if reactants.is_empty() {
        return Ok((reactants, mechanism));
    }

    let mut found_reactants: HashSet<String> = HashSet::new();
    let mut found_reactions: HashSet<String> = HashSet::new();
    const INERT_MARKERS: [&str; 7] = ["M", "M)", "(M", "(M)", "V", "(V", "V)"];

    loop {
        for (r_id, reaction_data) in db_object
            .iter()
            .filter(|(r_id, _)| !mechanism.contains(r_id.as_str()))
        {
            let reactants_db = reaction_data.get("reagents").ok_or_else(|| {
                KineticsError::InvalidReactionPayload {
                    reaction: r_id.clone(),
                    field: "reagents".into(),
                    reason: "relationship row is missing this field".into(),
                }
            })?;
            let products_db = reaction_data.get("products").ok_or_else(|| {
                KineticsError::InvalidReactionPayload {
                    reaction: r_id.clone(),
                    field: "products".into(),
                    reason: "relationship row is missing this field".into(),
                }
            })?;

            let reagents_ok = reactants_db
                .iter()
                .filter(|subs| !INERT_MARKERS.contains(&subs.as_str()))
                .all(|subs| reactants.contains(subs));

            if reagents_ok {
                found_reactants.extend(
                    products_db
                        .iter()
                        .filter(|product| !INERT_MARKERS.contains(&product.as_str()))
                        .cloned(),
                );
                found_reactions.insert(r_id.clone());
            }
        }

        if found_reactions.is_empty() {
            return Ok((reactants, mechanism));
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

    let (reactlibrary, reaction_db) =
        cached_catalog_pair().map_err(|error| error.with_context("loading mechanism catalogs"))?;
    let reactlibrary =
        reactlibrary
            .get(big_mech)
            .ok_or_else(|| KineticsError::MissingCatalogLibrary {
                catalog: "reaction parameters",
                library: big_mech.to_string(),
            })?;

    let reaction_db =
        reaction_db
            .get(big_mech)
            .ok_or_else(|| KineticsError::MissingCatalogLibrary {
                catalog: "reaction relationships",
                library: big_mech.to_string(),
            })?;
    let (mut reactants, mechanism) = parse_database(reaction_db, mechanism, reactants)
        .map_err(|error| error.with_context(format!("building mechanism `{big_mech}`")))?;

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
    let mechanism = vec_string_sorted(big_mech, &mechanism)?;

    let mut vec_of_reactions: Vec<String> = Vec::with_capacity(mechanism.len());
    let mut vec_of_reaction_value: Vec<Value> = Vec::with_capacity(mechanism.len());
    for react_num in &mechanism {
        let data_for_react = reactlibrary
            .get(react_num)
            .ok_or_else(|| KineticsError::MissingReaction(react_num.clone()))?;
        vec_of_reaction_value.push(data_for_react.to_owned());
        vec_of_reactions.push(reaction_equation(react_num, data_for_react)?.to_string());
    }
    Ok((
        mechanism,
        reactants,
        vec_of_reactions,
        vec_of_reaction_value,
    ))
}

#[cfg(test)]
mod tests {
    use super::{parse_database, reaction_equation, vec_string_sorted};
    use crate::Kinetics::error::KineticsError;
    use serde_json::json;
    use std::collections::{HashMap, HashSet};

    fn reaction(reagents: &[&str], products: &[&str]) -> HashMap<String, HashSet<String>> {
        HashMap::from([
            (
                "reagents".into(),
                reagents.iter().map(|item| (*item).into()).collect(),
            ),
            (
                "products".into(),
                products.iter().map(|item| (*item).into()).collect(),
            ),
        ])
    }

    #[test]
    fn graph_closure_is_forward_only_ignores_inert_markers_and_is_seeded() {
        let database = HashMap::from([
            ("1".into(), reaction(&["A", "(M)"], &["B", "M)"])),
            ("2".into(), reaction(&["B", "C"], &["D"])),
            ("3".into(), reaction(&["B"], &["C"])),
            ("4".into(), reaction(&["Z"], &["A"])),
        ]);

        let (species, selected) =
            parse_database(&database, HashSet::new(), HashSet::from(["A".into()])).unwrap();
        assert_eq!(
            species,
            HashSet::from(["A".into(), "B".into(), "C".into(), "D".into()])
        );
        assert_eq!(
            selected,
            HashSet::from(["1".into(), "2".into(), "3".into()])
        );

        let (empty_species, empty_selected) =
            parse_database(&database, HashSet::new(), HashSet::new()).unwrap();
        assert!(empty_species.is_empty());
        assert!(empty_selected.is_empty());
    }

    #[test]
    fn graph_closure_rejects_missing_relationship_fields_with_reaction_context() {
        let database = HashMap::from([("42".into(), HashMap::new())]);
        let error =
            parse_database(&database, HashSet::new(), HashSet::from(["A".into()])).unwrap_err();
        assert!(matches!(
            error,
            KineticsError::InvalidReactionPayload { reaction, field, .. }
                if reaction == "42" && field == "reagents"
        ));
    }

    #[test]
    fn numeric_id_ordering_preserves_original_strings_and_large_values() {
        let ids = vec![
            "10".into(),
            "0002".into(),
            "2".into(),
            "184467440737095516160".into(),
        ];
        assert_eq!(
            vec_string_sorted("fixture", &ids).unwrap(),
            vec!["0002", "2", "10", "184467440737095516160"]
        );
    }

    #[test]
    fn mechanism_id_ordering_rejects_non_numeric_ids_with_context() {
        let err = vec_string_sorted("fixture", &["R1".into()]).unwrap_err();
        assert!(matches!(
            err,
            crate::Kinetics::error::KineticsError::InvalidMechanismReactionId {
                mechanism,
                reaction_id
            } if mechanism == "fixture" && reaction_id == "R1"
        ));
    }

    #[test]
    fn equation_extraction_returns_string_content_not_json_encoding() {
        let payload = json!({"eq": "H2 + O2 => H2O"});
        assert_eq!(reaction_equation("R1", &payload).unwrap(), "H2 + O2 => H2O");
    }

    #[test]
    fn equation_extraction_rejects_non_string_payload_with_reaction_context() {
        let payload = json!({"eq": ["H2", "O2"]});
        assert!(matches!(
            reaction_equation("R7", &payload),
            Err(crate::Kinetics::error::KineticsError::InvalidReactionPayload {
                reaction,
                field,
                ..
            }) if reaction == "R7" && field == "eq"
        ));
    }
}
