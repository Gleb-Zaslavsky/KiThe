#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;
#[allow(non_snake_case)]
pub mod Examples;


use Examples::kinetics_examples::kin_examples;

pub fn main() {
    //
    let kintask: usize = 3;
    kin_examples(kintask);

    
}
