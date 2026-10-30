/// 0.1.1
/// ru
///  в проекте используются шорткаты для названий реакций в стиле
/// "сокращенное название библиотеки реакции"+"номер реакции в библиотеке",
///  например C_10 означает реакция '10' в библиотеке CANTERA.
/// То же означает C10 и Сantera_10
/// eng
/// The project uses shortcuts for reaction names in the style
/// of "short name of the reaction library" + "reaction number in the library"
/// for example C_10 means reaction '10' in the CANTERA library.
/// The same means C10 and Сantera_10
use std::collections::HashMap;

use crate::Kinetics::error::{KineticsError, KineticsResult};

/// Parse `<library alias><optional underscore><integer reaction id>`.
pub fn decipher_shortcut_name(short: &str) -> KineticsResult<(String, String)> {
    let input = short.trim();
    let digit_start = input
        .find(|character: char| character.is_ascii_digit())
        .ok_or_else(|| invalid_shortcut(short, "expected a reaction number"))?;
    let (mut prefix, reaction_id) = input.split_at(digit_start);
    if prefix.ends_with('_') {
        prefix = &prefix[..prefix.len() - 1];
    }
    if prefix.is_empty()
        || !prefix
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
        || reaction_id.is_empty()
        || !reaction_id.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(invalid_shortcut(
            short,
            "expected an ASCII library alias followed by an integer reaction id",
        ));
    }

    let normalized = prefix.to_ascii_lowercase().replace('_', "");
    let library = match normalized.as_str() {
        "c" | "ca" | "cantera" => "Cantera".to_string(),
        "n" | "nuig" => "NUIG".to_string(),
        "a" | "aramco" => "Aramco".to_string(),
        "b" | "beckstead" => "Beckstead".to_string(),
        "bc" | "becksteadc" => "Beckstead_c".to_string(),
        _ => prefix.to_string(),
    };
    Ok((library, reaction_id.to_string()))
}

fn invalid_shortcut(input: &str, reason: &str) -> KineticsError {
    KineticsError::InvalidReactionShortcut {
        input: input.to_string(),
        reason: reason.to_string(),
    }
}

pub fn decipher_vector_of_shortcuts(
    vec_of_shortcuts: Vec<&str>,
) -> KineticsResult<HashMap<String, Vec<String>>> {
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    for shortcut in vec_of_shortcuts {
        let (library_name, reaction_number) = decipher_shortcut_name(shortcut)?;
        map.entry(library_name).or_default().push(reaction_number);
    }
    Ok(map)
}

pub fn decipher_vector_of_shortcuts_to_pairs(
    vec_of_shortcuts: Vec<&str>,
) -> KineticsResult<Vec<(String, String)>> {
    let mut vec_of_pairs = Vec::new();
    for shortcut in vec_of_shortcuts {
        let (library_name, reaction_number) = decipher_shortcut_name(shortcut)?;
        vec_of_pairs.push((library_name, reaction_number));
    }
    Ok(vec_of_pairs)
}
pub fn cipher_vector_of_shortcuts(map: HashMap<String, String>) -> Vec<String> {
    let mut shortcuts: Vec<_> = map
        .into_iter()
        .map(|(library, reaction)| format!("{library}_{reaction}"))
        .collect();
    shortcuts.sort();
    shortcuts
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashMap, HashSet};
    #[test]
    fn test_decipher_shortcut_name() {
        assert_eq!(
            decipher_shortcut_name("C_10").unwrap(),
            ("Cantera".to_string(), "10".to_string())
        );
        assert_eq!(
            decipher_shortcut_name("C10").unwrap(),
            ("Cantera".to_string(), "10".to_string())
        );
        assert_eq!(
            decipher_shortcut_name("Cantera_10").unwrap(),
            ("Cantera".to_string(), "10".to_string())
        );

        assert_eq!(
            decipher_shortcut_name("A10").unwrap(),
            ("Aramco".to_string(), "10".to_string())
        );

        assert_eq!(
            decipher_shortcut_name("N_10").unwrap(),
            ("NUIG".to_string(), "10".to_string())
        );

        assert_eq!(
            decipher_shortcut_name("Bc_10").unwrap(),
            ("Beckstead_c".to_string(), "10".to_string())
        );
        assert_eq!(
            decipher_shortcut_name("B10").unwrap(),
            ("Beckstead".to_string(), "10".to_string())
        );
        assert_eq!(
            decipher_shortcut_name("NoSuchLib_10").unwrap(),
            ("NoSuchLib".to_string(), "10".to_string())
        );
        assert!(matches!(
            decipher_shortcut_name("junk-C10"),
            Err(KineticsError::InvalidReactionShortcut { .. })
        ));
        assert!(matches!(
            decipher_shortcut_name("C_10tail"),
            Err(KineticsError::InvalidReactionShortcut { .. })
        ));
        assert!(matches!(
            decipher_shortcut_name("invalid"),
            Err(KineticsError::InvalidReactionShortcut { .. })
        ));
    }

    #[test]
    fn test_decipher_vector_of_shortcuts() {
        let vec_of_shortcuts = vec!["C_10", "A10", "n_10"];
        let expected_map = HashMap::from([
            ("Cantera".to_string(), vec!["10".to_string()]),
            ("Aramco".to_string(), vec!["10".to_string()]),
            ("NUIG".to_string(), vec!["10".to_string()]),
        ]);

        assert_eq!(
            decipher_vector_of_shortcuts(vec_of_shortcuts).unwrap(),
            expected_map
        );
    }
    #[test]
    fn test_cipher_vector_of_shortcuts() {
        let map = vec![
            ("Cantera".to_string(), "10".to_string()),
            ("ARAMCO".to_string(), "10".to_string()),
            ("NUIG".to_string(), "10".to_string()),
        ]
        .into_iter()
        .collect::<HashMap<String, String>>();

        let expected_vec = vec![
            "Cantera_10".to_string(),
            "ARAMCO_10".to_string(),
            "NUIG_10".to_string(),
        ];
        let res: HashSet<String> = cipher_vector_of_shortcuts(map).into_iter().collect();
        let exp: HashSet<String> = expected_vec.into_iter().collect();
        assert_eq!(res, exp);
    }
}
