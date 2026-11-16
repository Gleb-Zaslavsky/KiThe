//! Numerical methods used by experimental-kinetics datasets.
//!
//! This module groups fitting, smoothing/resampling filters, and statistical
//! quality measures. Dataset and series orchestration remain in the parent module.

pub mod exp_kinetics_smooth_filter;
pub mod fitting;
pub mod ndarray_statistics;

#[cfg(test)]
mod filters_smoothing_test;
