use std::collections::HashMap;
use std::f64::consts::E;
const R: f64 = 1.987; // кал/(K·моль).
fn Cp(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    R*(a + b * t + c * t.powi(2) + d * t.powi(3) + e * t.powi(4) )
   
}
fn dh(t: f64, a: f64, b: f64, c: f64, d: f64, e: f64, f: f64) -> f64 {
    R * t * (a + b * t / 2.0 + c * t.powi(2) / 3.0 + d * t.powi(3) / 4.0 + e * t.powi(4) / 5.0 + f / t)
}
fn ds(t : f64, a: f64, b: f64, c: f64, d: f64, e: f64, g: f64) -> f64 {
    R * (a * t.ln() + b * t + c * t.powi(2) / 2.0 + d * t.powi(3) / 3.0 + e * t.powi(4) / 4.0 + g);
}

// src\DBhandlers\NASAdata.rs
pub struct NASAdata {
    pub data_substance:  Vec<f64>,
    pub Cp: f64,
    pub dh: f64,
    pub ds: f64,
}


impl NASAdata {
    pub fn new() -> Self {  Self { data_substance: HashMap::new() } }

    pub fn set_data_substance(&mut self, data_substance: HashMap<String, Vec<f64>>) {
        self.data_substance = data_substance;
    }

    pub fn calculate_Cp_dH_dS(&mut self, t: f64)  {
    
    let c_data = self.data_substance;

 

    if c_data.len() == 25 {
        let t1 = c_data[0];
        let t2 = c_data[1];
        let t3 = c_data[2];
        let t4 = c_data[3];
        if t1 <= t && t <= t2 {
            let a = c_data[4];
            let b = c_data[5];
            let c = c_data[6];
            let d = c_data[7];
            let e = c_data[8];
            let f = c_data[9];
            let g = c_data[10];

        } else if t2 < t && t <= t3 {
            let a = c_data[11];
            let b = c_data[12];
            let c = c_data[13];
            let d = c_data[14];
            let e = c_data[15];
            let f = c_data[16];
            let g = c_data[17];

        } else if t3 < t && t < t4 {
            let a = c_data[18];
            let b = c_data[19];
            let c = c_data[20];
            let d = c_data[21];
            let e = c_data[22];
            let f = c_data[23];
            let g = c_data[24];

        }//
    }//  end of if c_data.len() == 25
    else if c_data.len() == 17 {
        let t1 = c_data[0];
        let t2 = c_data[1];
        let t3 = c_data[2];
    
        if t1 <= t && t <= t2 {
            let (a, b, c, d, e, f, g) = (c_data[3], c_data[4], c_data[5], c_data[6], c_data[7], c_data[8], c_data[9]);
        } else if t2 < t && t <= t3 {
            let (a, b, c, d, e, f, g) = (c_data[10], c_data[11], c_data[12], c_data[13], c_data[14], c_data[15], c_data[16]);
        }  } 
    else if c_data.len() == 9 {
        let t1 = c_data[0];
        let t2 = c_data[1];
    
        if t1 <= t && t <= t2 {
            let (a, b, c, d, e, f, g) = (c_data[2], c_data[3], c_data[4], c_data[5], c_data[6], c_data[7], c_data[8]);
        } } 

    self. Cp = Cp(t, a, b, c, d, e) ;
     self.dh = dh(t, a, b, c, d, e, f);
    self.ds = ds(t, a, b, c, d, e, g);

}
}

