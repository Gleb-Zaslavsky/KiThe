pub mod column_provenance;
pub mod exp_engine_api;
pub mod exp_kinetics_column_manipulation;
#[path = "experimental_kinetics/experiment_series/experiment_series2.rs"]
pub mod experiment_series2;
#[path = "experimental_kinetics/experiment_series/experiment_series3.rs"]
pub mod experiment_series3;
#[path = "experimental_kinetics/experiment_series/experiment_series_main.rs"]
pub mod experiment_series_main;
#[cfg(test)]
#[path = "experimental_kinetics/experiment_series/experiment_series_main_test.rs"]
pub(crate) mod experiment_series_main_test;
pub mod math;
pub mod one_experiment_dataset;
#[cfg(test)]
pub(crate) mod one_experiment_dataset_test;

pub mod testing_mod;
//pub mod LSQSplines2;

#[cfg(test)]
mod invariant_tests;
pub mod kinetic_methods;
pub mod kinetic_methods2;
#[cfg(test)]
pub(crate) mod kinetic_methods_tests;
