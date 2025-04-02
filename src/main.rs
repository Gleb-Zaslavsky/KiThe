#[allow(non_snake_case)]
pub mod Examples;
#[allow(non_snake_case)]
pub mod Kinetics;
#[allow(non_snake_case)]
pub mod ReactorsBVP;
#[allow(non_snake_case)]
pub mod Thermodynamics;
#[allow(non_snake_case)]
pub mod Utils;
#[allow(unused_imports)]
use Examples::NIST_examples::NIST_examples;
#[allow(unused_imports)]
use Examples::SubstanceDataCollecting::collecting_thermo_data;
#[allow(unused_imports)]
use Examples::kinetics_examples::kin_examples;
#[allow(unused_imports)]
use Examples::thermo_examples::thermo_examples;

pub fn main() {
    //
    #[allow(unused_variables)]
    let task: usize = 2;
    //kin_examples(3);
    //  thermo_examples(task);
    // NIST_examples(task);
    collecting_thermo_data(task);
}
