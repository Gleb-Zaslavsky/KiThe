//v 0.1.1
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
//use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::time::Instant;
//use json_to_table::json_to_table;
#[warn(unused_imports)]
#[warn(unused_variables)]
#[warn(non_snake_case)]

fn convert(
    input: &HashMap<String, HashMap<String, HashMap<String, HashSet<String>>>>,
) -> HashMap<&str, HashMap<&str, HashMap<&str, HashSet<&str>>>> {
    let mut output: HashMap<&str, HashMap<&str, HashMap<&str, HashSet<&str>>>> = HashMap::new();

    for (k1, v1) in input.iter() {
        let mut nested_map1: HashMap<&str, HashMap<&str, HashSet<&str>>> = HashMap::new();

        for (k2, v2) in v1.iter() {
            let mut nested_map2: HashMap<&str, HashSet<&str>> = HashMap::new();

            for (k3, v3) in v2.iter() {
                let mut set: HashSet<&str> = HashSet::new();

                for value in v3.iter() {
                    set.insert(value);
                }

                nested_map2.insert(k3, set);
            }

            nested_map1.insert(k2, nested_map2);
        }

        output.insert(k1, nested_map1);
    }

    output
}
pub fn hashset<'a>(data: &'a Vec<&str>) -> HashSet<&'a str> {
    HashSet::from_iter(data.iter().cloned())
}

fn vec_string_sorted<'a>(vec_strings: &'a Vec<String>) -> Vec<String> {
    let mut vec_ints: Vec<i32> = vec_strings
        .iter()
        .map(|s| s.parse::<i32>().unwrap())
        .collect();
    vec_ints.sort();
    let vec_sorted_strings: Vec<String> = vec_ints.iter().map(|i| i.to_string()).collect();
    vec_sorted_strings
}
fn convert_vec_string_to_vec_str(vec: &Vec<String>) -> Vec<&str> {
    vec.iter().map(|s| s.as_str()).collect()
}
pub fn parse_database<'a>(
    db_object_: &'a HashMap<&str, HashMap<&str, HashSet<&'a str>>>,
    mut mechanism: HashSet<&'a str>,
    mut reactants: HashSet<&'a str>,
) -> (HashSet<String>, HashSet<String>) {
    let mut found_reactants: HashSet<&str> = HashSet::new();
    let mut found_reactions: HashSet<&str> = HashSet::new();
    let db_object = &db_object_;
    let db_object_mut = &mut db_object_.clone();
    let keys_to_replace = HashSet::from(["M", "M)", "(M", "(M)"]);
    loop {
        let all_reactions = db_object.keys();
        let all_reactions_ = all_reactions.filter(|r_id| !mechanism.contains(*r_id));
        for r_id in all_reactions_ {
            //делаем изменяемым hashset с веществами, чтобы можно было выбросить оттуда keys_to_replace
            let reactants_db = &mut db_object_mut
                .get_mut(r_id)
                .expect("REASON")
                .get_mut("reagents")
                .expect("REASON");
            // сохраняем все значения, которые не входят в keys_to_replace
            reactants_db.retain(|subs| !keys_to_replace.contains(subs));
            if reactants_db.is_subset(&reactants) {
                let products_db = db_object[r_id]["products"].iter();
                found_reactants.extend(products_db);
                found_reactions.insert(r_id);
            }
        }
        if found_reactions.is_empty() {
            return (
                reactants.iter().cloned().map(String::from).collect(),
                mechanism.iter().cloned().map(String::from).collect(),
            );
        }
        mechanism.extend(&found_reactions);
        reactants.extend(&found_reactants);
        found_reactants.clear();
        found_reactions.clear();
    }
}
pub fn mechfinder(
    big_mech: &str,
    vec: Vec<&str>,
) -> (Vec<String>, Vec<String>, Vec<String>, Vec<Value>) {
    // v4: задача "NUIG"/ ("O", "NH3", "NO") всего - 90 мс  парсинг данных 37 мс (в т.ч.convert_hashmap - 8, разбор  reaction_db - 8) , расчет 45, вывод 7
    let before = Instant::now(); //начало отсчета времени выполнения
    let search_s = hashset(&vec);
    let mech: HashSet<&str> = HashSet::new();
    let reactants: HashSet<&str> = search_s.clone();
    //
    let mut file = File::open("Reactbase.json").unwrap();
    let mut file_contents = String::new();
    file.read_to_string(&mut file_contents).unwrap();
    // библиотека:{ номер реакции :{ данные реакции }}
    let reactlibrary: HashMap<String, HashMap<String, Value>> =
        serde_json::from_str::<HashMap<String, HashMap<String, _>>>(&file_contents).unwrap();
    let reactlibrary = &reactlibrary[big_mech];

    //
    let mut file = File::open("dict_reaction.json").unwrap();
    let mut file_contents = String::new();
    file.read_to_string(&mut file_contents).unwrap();
    // библиотека:{ номер реакции :{ "reagents"/"products":{}: }}
    let reaction_db: HashMap<String, HashMap<String, HashMap<String, HashSet<String>>>> =
        serde_json::from_str(&file_contents).unwrap();

    let reaction_db: HashMap<&str, HashMap<&str, HashMap<&str, HashSet<&str>>>> =
        convert(&reaction_db);

    let (mut reactants, mechanism) = parse_database(&reaction_db[big_mech], mech, reactants);
    let keys_to_replace = HashSet::from(["M", "M)", "(M", "(M)"]);
    reactants.retain(|subs| !keys_to_replace.contains(subs.as_str()));
    let mut reactants = reactants.into_iter().collect::<Vec<_>>();
    reactants.sort();
    let mechanism: Vec<_> = mechanism.into_iter().collect();
    let mechanism_: Vec<String> = vec_string_sorted(&mechanism);
    let mechanism_: Vec<&str> = convert_vec_string_to_vec_str(&mechanism_);
    //  println!("ВЕЩЕСТВА: {:?}", &reactants);
    //   println!("РЕАКЦИИ:: {:?}",  &mechanism_);
    println!("найдено веществ, {}", &reactants.len());
    println!("найдено реакций, {}", &mechanism_.len());
    let mut vec_of_reactions: Vec<String> = Vec::new();
    let mut vec_of_reaction_value: Vec<Value> = Vec::new();
    for react_num in &mechanism {
        let data_for_react = &reactlibrary[react_num];
        vec_of_reaction_value.push(data_for_react.to_owned());
        //  let table = json_to_table(data_for_react).to_string();
        let json_string = serde_json::to_string(data_for_react).unwrap();
        vec_of_reactions.push(json_string);
        //  println!("{}", json_string);
        //  println!("таблица, {}", table)
    }
    let time_elapsed = before.elapsed();
    println!("Elapsed time: {:.2?}", time_elapsed);
    return (
        mechanism,
        reactants,
        vec_of_reactions,
        vec_of_reaction_value,
    );
}
