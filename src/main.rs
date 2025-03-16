#[allow(non_snake_case)]
pub mod Examples;
#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;

use Examples::kinetics_examples::kin_examples;
use Examples::thermo_examples::thermo_examples;
use Examples::NIST_examples::NIST_examples;

pub fn main() {
    //
    let task: usize = 0;
    //thermo_examples(task);
    NIST_examples(task);
}
