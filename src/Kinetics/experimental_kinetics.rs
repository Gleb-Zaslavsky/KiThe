pub mod column_provenance;
pub mod exp_engine_api;
pub mod exp_kinetics_column_manipulation;
pub mod exp_kinetics_smooth_filter;
#[path = "experimental_kinetics/experiment_series/experiment_series2.rs"]
pub mod experiment_series2;
#[path = "experimental_kinetics/experiment_series/experiment_series3.rs"]
pub mod experiment_series3;
#[path = "experimental_kinetics/experiment_series/experiment_series_main.rs"]
pub mod experiment_series_main;
#[path = "experimental_kinetics/experiment_series/experiment_series_main_test.rs"]
pub mod experiment_series_main_test;
mod filters_smoothing_test;
pub mod ndarray_statistics;
pub mod one_experiment_dataset;
pub mod one_experiment_dataset_test;

pub mod testing_mod;
//pub mod LSQSplines2;

pub mod fitting;
pub mod kinetic_methods;
pub mod kinetic_methods2;
pub mod kinetic_methods_tests;
