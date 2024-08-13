use std::collections::{HashMap, HashSet};
use serde::{Deserialize, Serialize};
use std::f32;

const R:f32 = 8.314;
#[derive(Debug, Deserialize, Serialize)]
pub struct  ElementaryStruct {
    pub  r#type:String,
    pub Arrenius: Vec<f32>,
    pub eq: String,
   // react:HashMap<String, i32>,
    pub react: Option<  HashMap<String, i32>, >,
    pub reactnumber: Option<  i32 >,
 
}
impl ElementaryStruct {
    pub fn new(    r#type:String, Arrenius: Vec<f32>, eq: String,react: Option<  HashMap<String, i32>, >, reactnumber: Option<  i32 >,) -> Self {
        Self {
            r#type,
            Arrenius,
            eq,
            react,
            reactnumber
        }
    }
    pub fn K_const(&self, Temp: f32 ) -> f32 {
        let A:&f32 = &self.Arrenius[0];
        let n:&f32 = &self.Arrenius[1];
        let E:&f32 = &self.Arrenius[2];
        let K_const_: f32 =  A*Temp.powf(*n)*f32::exp(-E/(Temp*R));
        return K_const_
    }
}
#[derive(Debug, Deserialize, Serialize)]
pub struct FalloffStruct {
    pub  r#type:String,
    pub Arrenius: Vec<f32>,
    pub eq: String,
}

pub struct PressureStruct {
    
}

pub struct ThreeBodyStruct {
    
}