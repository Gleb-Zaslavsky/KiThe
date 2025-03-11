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

pub fn main() {
    //
    let task: usize = 2;
    thermo_examples(task);
}
