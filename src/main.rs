#[allow(non_snake_case)]
pub mod Examples;
#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;
#[allow(unused_imports)]
use Examples::NIST_examples::NIST_examples;
#[allow(unused_imports)]
use Examples::kinetics_examples::kin_examples;
#[allow(unused_imports)]
use Examples::thermo_examples::thermo_examples;

pub fn main() {
    //
    let task: usize = 5;
    thermo_examples(task);
    // NIST_examples(task);
}
