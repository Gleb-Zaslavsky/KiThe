use super::SimpleReactorBVP::{FastElemReact, R_G, SimpleReactorTask};
use crate::ReactorsBVP::reactor_BVP_utils::ScalingConfig;
use std::collections::HashMap;
#[allow(dead_code)]
/// Build a compact HMX-style reactor fixture that exercises the full BVP pipeline.
pub(crate) fn compact_hmx_reactor() -> SimpleReactorTask {
    let eq = "HMX=>HMXprod".to_string();
    let c_p = 0.35 * 4.184 * 1000.0;
    let lambda_eff = 0.07;
    let n = 0.0;
    let m = 0.077 * (1e6_f64 / 1e5_f64).powf(0.748) / 1e2;
    let a = 1.3e5;
    let e = 5000.0 * 4.184;
    let t0 = 600.0;
    let t_scale = 600.0;
    let p = 1e6;
    let tm = 1500.0;
    let hmx = HashMap::from([
        ("H".to_string(), 4),
        ("N".to_string(), 8),
        ("C".to_string(), 8),
        ("O".to_string(), 8),
    ]);
    let hmxprod = HashMap::from([
        ("H".to_string(), 6),
        ("C".to_string(), 1),
        ("O".to_string(), 1),
    ]);
    let groups = Some(HashMap::from([
        ("HMX".to_string(), hmx.clone()),
        ("HMXprod".to_string(), hmxprod),
    ]));

    let mut reactor = SimpleReactorTask::new();
    let reactions = vec![FastElemReact {
        eq,
        A: a,
        n,
        E: e,
        Q: 3000.0 * 1e3 / 100.0,
    }];

    reactor
        .fast_react_set(reactions)
        .expect("valid elementary reaction should be accepted");
    reactor.kindata.substances = vec!["HMX".to_string(), "HMXprod".to_string()];
    reactor.kindata.groups = groups;

    let diffusion = {
        let ro0 = 34.2e-3 * p / (R_G * t0);
        let d = lambda_eff / (c_p * ro0);
        HashMap::from([("HMX".to_string(), d), ("HMXprod".to_string(), d)])
    };
    let boundary_condition = HashMap::from([
        ("HMX".to_string(), 1.0 - 1e-3),
        ("HMXprod".to_string(), 1e-3),
        ("T".to_string(), t0),
    ]);
    let thermal_effects = vec![3000.0 * 1e3 / 100.0];
    let scaling = ScalingConfig::new(t_scale, 9e-4, t_scale);
    reactor.set_parameters(
        thermal_effects,
        p,
        tm,
        c_p,
        boundary_condition,
        lambda_eff,
        diffusion,
        m,
        scaling,
    );
    reactor.M = 34.2 / 1000.0;
    reactor
}
