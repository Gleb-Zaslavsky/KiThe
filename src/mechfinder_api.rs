mod mechfinder;
pub mod kinetics; //pub!
use kinetics::{ ElementaryStruct, FalloffStruct, PressureStruct, ThreeBodyStruct};

use serde_json::{Value, Map, Number};
use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};
use std::f32;



#[warn(unused_imports)]
enum ReactionTypes {
    elem,
    falloff,
    pressure,
    threebody
}

/* Нам необходимо расширить существующую струкутуру, полученную парсингом json чтобы добавить */

fn parse_kinetic_data( big_mech:&str,  vec_of_reactions:&Vec<String>,  vec_of_reaction_value:Vec<Value>)->() {
    let mut ReactionDataHash:   Map<String, Value>  = Map::new();
    //let mut parsed_reactions: Vec<_> = Vec::new();
    for reaction_record in &vec_of_reaction_value {
        if let Some(j) = &vec_of_reaction_value.iter().position(|x| x == reaction_record) {
            
        let react_id:String = vec_of_reactions[*j].to_owned();
        println!("{:?}", &reaction_record);
        let react_type_value:&Value = &reaction_record["type"];
       // println!("{:?}", &react_type_value);
        let react_type_value:String =  serde_json::from_value(react_type_value.to_owned()).unwrap();
       // print!("тип реакции {:?}", &react_type_value);
        match react_type_value.as_str() {

            "elem" =>{
                 // парсим json запись в структуру, конструкция if let Ok нужна для обработки ошибки при записи в структуру, Т.К. 
                 // в структуре есть поля которых нет в json
                 if let Ok(mut elemstruct) = serde_json::from_value::<ElementaryStruct>(reaction_record.clone()) {
                    elemstruct.react = Some(HashMap::new());
                    elemstruct.reactnumber = Some(0);
                     println!("{:?}", &elemstruct);   
                     let v = serde_json::to_value(&elemstruct).unwrap();
                     ReactionDataHash.insert(react_id.clone(), v );
                     /* 
                     get value syntax
                     let v:ElementaryStruct = serde_json::from_value::<ElementaryStruct>(ReactionDataHash[&react_id].clone()).unwrap();
                     print!("value {:?}", v);
                     println!("K_const {:?}", v.K_const(298.15) );
                    */
                } 
                     }
            "falloff" =>{
                if let Ok(mut faloffstruct) = serde_json::from_value::<FalloffStruct>(reaction_record.clone()) {
                  //  elemstruct.react = Some(HashMap::new());
                   // elemstruct.reactnumber = Some(0);
                     println!("{:?}", &faloffstruct);   
                     let v = serde_json::to_value(&faloffstruct).unwrap();
                     ReactionDataHash.insert(react_id, v );
                   // parsed_reactions.push(elemstruct);
                    
                } 

             }
            "pressure" =>{

             }
            "threebody" =>{

             }
            _ => {  print!("")}
    }; // match

    } // if let Some
  //  ReactionDataHash.insert(react_id, elemstruct);
} // for
}
#[derive(Debug)]
 pub struct Mechanism_search {
    pub task_substances:Vec<String>, 
    pub task_library:String, 
    pub mechanism: Vec<String>, 
    pub reactants: Vec<String>,
    pub vec_of_reactions: Vec<String>, 
   

}
 // 
impl Mechanism_search {
    pub fn new(task_substances:Vec<String>, task_library:String, mechanism:Vec<String>, reactants:Vec<String>, vec_of_reactions:Vec<String>)-> Self { 
        Self{
            task_substances,
            task_library,
            mechanism,
            reactants,
            vec_of_reactions,
        } 
        }
    
    pub fn default() -> Self {    
        Self{
            task_substances:Vec::new(), 
            task_library: String::new(),
            mechanism: Vec::new(),
            reactants: Vec::new(),
            vec_of_reactions: Vec::new(),
        }
    }
    pub fn mechfinder_api(&mut self)-> ( Vec<String>, Vec<String>, Vec<String>,) {
        /* 
        let tuple = [ "O", "NH3", "NO", "O2", "N2", "N2O", "CO", "C"];
        O,NH3,NO,O2,N2,N2O,CO,C
        let big_mech = "NUIG".to_string();
            let vec: Vec<&str> =  tuple.into_iter().collect();
        */
    

        let vec: Vec<&str> = self.task_substances.iter().map(|s| s.as_str()).collect();
        let big_mech = self.task_library.clone();
        println!("задание {:?}, библиотека {:?}", &big_mech, &vec);
        let (mechanism, reactants, vec_of_reactions,  vec_of_reaction_value) 
        = mechfinder::mechfinder(&big_mech, vec);
       // println!("mechanism {:?}", &mechanism);
         parse_kinetic_data(&big_mech,  &vec_of_reactions, vec_of_reaction_value.clone());  // парсим данные о реакциях
        self.mechanism = mechanism;
        self.reactants = reactants;
        self.vec_of_reactions = vec_of_reactions;
        return (self.mechanism.to_owned(), self.reactants.to_owned(), self.vec_of_reactions.to_owned()) 
        
    }
}
