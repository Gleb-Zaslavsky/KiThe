use super::*;
use crate::ReactorsIVP::SimpleReactorIVP::SimpleReactorTask;
use RustedSciThe::command_interpreter::task_parser::{DocumentMap, DocumentParser};
use RustedSciThe::command_interpreter::task_parser_ivp::{IvpMethodSpec, Lsode2TaskExecutionSpec};

fn parse_document_for_ivp(input: &str) -> DocumentMap {
    let mut parser = DocumentParser::new(input.to_string());
    parser.parse_document().expect("IVP document should parse");
    parser
        .get_result()
        .expect("IVP document map should exist")
        .clone()
}

#[test]
fn test_parse_reactor_ivp_document_delegates_solver_settings() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
problem_name: CondensedBurn
problem_description: Condensed-phase test document
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

scaling
dT: 100.0
T_scale: 75.0
L: 0.02

initial_conditions
T: 900.0
q: 0.0
A: 1.0

solver_options
lsode2_method_family: auto
lsode2_symbolic_assembly: AtomView
lsode2_symbolic_execution: AOT
lsode2_aot_toolchain: c_tcc
lsode2_aot_profile: release
lsode2_linear_structure: sparse
lsode2_native_execution: faithful_bdf_solve
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("reactor IVP document should parse successfully");

    assert_eq!(spec.problem_name.as_deref(), Some("CondensedBurn"));
    assert_eq!(
        spec.problem_description.as_deref(),
        Some("Condensed-phase test document")
    );
    assert_eq!(spec.physical.ro, 1800.0);
    assert_eq!(spec.physical.Cp, 1200.0);
    assert_eq!(spec.physical.Lambda, 0.35);
    assert_eq!(spec.physical.m, 0.02);
    assert_eq!(spec.physical.L, 0.02);
    assert_eq!(spec.scaling.L, 0.02);
    assert_eq!(spec.scaling.dT, 100.0);
    assert_eq!(spec.scaling.T_scale, 75.0);
    assert_eq!(spec.initial_conditions["T"], 900.0);
    assert_eq!(spec.initial_conditions["q"], 0.0);
    assert_eq!(spec.initial_conditions["A"], 1.0);
    assert_eq!(spec.thermal_effects, vec![125000.0]);
    assert!(!spec.stop_condition.enabled);
    assert!(spec.stop_condition.species.is_none());
    assert!(spec.stop_condition.threshold.is_none());
    assert_eq!(spec.solver_settings.solver.method, IvpMethodSpec::Lsode2);

    let lsode2 = spec
        .solver_settings
        .solver_options
        .lsode2
        .as_ref()
        .expect("LSODE2-specific options should be present");
    assert!(matches!(
        lsode2.symbolic_execution,
        Some(Lsode2TaskExecutionSpec::Aot { .. })
    ));
}

#[test]
fn test_canonical_reactor_ivp_template_is_parser_ready() {
    let template = build_canonical_condensed_reactor_ivp_task_template();
    let spec = parse_reactor_ivp_document_from_str(&template)
        .expect("canonical reactor IVP template should parse");

    assert_eq!(spec.problem_name.as_deref(), Some("CondensedBurn"));
    assert_eq!(
        spec.problem_description.as_deref(),
        Some("Canonical condensed-phase IVP template")
    );
    assert_eq!(spec.physical.ro, 1800.0);
    assert_eq!(spec.physical.Cp, 1200.0);
    assert_eq!(spec.physical.Lambda, 0.35);
    assert_eq!(spec.physical.m, 0.02);
    assert_eq!(spec.physical.L, 0.02);
    assert_eq!(spec.scaling.L, 0.02);
    assert_eq!(spec.initial_conditions["T"], 900.0);
    assert_eq!(spec.initial_conditions["q"], 0.0);
    assert_eq!(spec.initial_conditions["A"], 1.0);
    assert_eq!(spec.initial_conditions["B"], 0.0);
    assert_eq!(spec.thermal_effects, vec![125000.0]);
    assert_eq!(spec.solver_settings.solver.method, IvpMethodSpec::Lsode2);
    assert!(spec.solver_settings.solver_options.lsode2.is_some());
}

#[test]
fn test_parse_reactor_ivp_document_delegates_integer_and_optional_solver_settings() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0

solver_options
max_iterations: 42
rtol: 1e-6
atol: 1e-8
max_step: 1e-3
first_step: Some(1e-6)
parallel: false
lsode2_method_family: auto
lsode2_symbolic_assembly: AtomView
lsode2_symbolic_execution: LambdifyExpr
lsode2_linear_structure: sparse
"#;

    let spec =
        parse_reactor_ivp_document_from_str(input).expect("typed solver settings should parse");

    assert_eq!(spec.solver_settings.solver_options.max_iterations, Some(42));
    assert_eq!(spec.solver_settings.solver_options.rtol, Some(1e-6));
    assert_eq!(spec.solver_settings.solver_options.atol, Some(1e-8));
    assert_eq!(spec.solver_settings.solver_options.max_step, Some(1e-3));
    assert_eq!(spec.solver_settings.solver_options.first_step, Some(1e-6));
    assert_eq!(spec.solver_settings.solver_options.parallel, Some(false));
    assert!(spec.solver_settings.solver_options.lsode2.is_some());
}

#[test]
fn test_parse_reactor_ivp_document_reads_stop_condition_section() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0
B: 0.0

stop_condition
enabled: true
species: A
threshold: 1e-4
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("stop-condition document should parse successfully");

    assert!(spec.stop_condition.enabled);
    assert_eq!(spec.stop_condition.species.as_deref(), Some("A"));
    assert_eq!(spec.stop_condition.threshold, Some(1.0e-4));
}

#[test]
fn test_parse_reactor_ivp_document_requires_physics_section() {
    let input = r#"
task
solver: IVP
method: LSODE2

initial_conditions
T: 900.0
q: 0.0
A: 1.0
"#;

    let error = parse_reactor_ivp_document_from_str(input)
        .expect_err("missing physics section should be rejected");
    assert!(error.to_string().contains("process_conditions"));
}

#[test]
fn test_apply_reactor_ivp_document_copies_physics_state() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("reactor IVP document should parse successfully");
    let mut task = SimpleReactorTask::new();

    spec.apply_to_task(&mut task)
        .expect("physics contract should apply to the task");

    assert_eq!(task.ro, 1800.0);
    assert_eq!(task.Cp, 1200.0);
    assert_eq!(task.Lambda, 0.35);
    assert_eq!(task.m, 0.02);
    assert_eq!(task.L, 0.02);
    assert_eq!(task.thermal_effects, vec![125000.0]);
    assert_eq!(task.initial_conditions["T"], 900.0);
    assert_eq!(task.initial_conditions["q"], 0.0);
    assert_eq!(task.initial_conditions["A"], 1.0);
}

#[test]
fn test_validate_against_species_rejects_missing_initial_entries() {
    let document = parse_document_for_ivp(
        r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0
"#,
    );

    let spec = parse_reactor_ivp_document(&document).expect("document should parse");
    let error = spec
        .validate_against_species(&["B".to_string()])
        .expect_err("missing species should be rejected");
    assert!(error.to_string().contains("B"));
}

#[test]
fn test_validate_against_species_rejects_unknown_stop_condition_species() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0

stop_condition
enabled: true
species: UnknownFuel
threshold: 1e-4
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("stop-condition document should parse successfully");

    let error = spec
        .validate_against_species(&["A".to_string()])
        .expect_err("unknown stop-condition species should be rejected");
    assert!(error.to_string().contains("UnknownFuel"));
}

#[test]
fn test_validate_document_is_pure_and_does_not_mutate_task() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

initial_conditions
T: 900.0
q: 0.0
A: 1.0
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("reactor IVP document should parse successfully");
    let mut task = SimpleReactorTask::new();
    task.set_density(1550.0)
        .expect("density should be accepted");
    task.m = 0.03;
    task.L = 0.07;

    let snapshot_before = (task.ro, task.m, task.L, task.initial_conditions.clone());
    let result = spec.validate(&["A".to_string(), "B".to_string()]);

    assert!(result.is_err(), "missing B should make validation fail");
    assert_eq!(task.ro, snapshot_before.0);
    assert_eq!(task.m, snapshot_before.1);
    assert_eq!(task.L, snapshot_before.2);
    assert_eq!(task.initial_conditions, snapshot_before.3);
}

#[test]
fn test_reactor_ivp_document_roundtrip_preserves_solver_and_physics_contract() {
    let input = r#"
task
solver: IVP
method: LSODE2

process_conditions
problem_name: RoundTripCondensed
problem_description: Roundtrip document
ro: 1800.0
Cp: 1200.0
Lambda: 0.35
m: 0.02
L: 0.02
thermal_effects: [125000.0]

scaling
dT: 100.0
T_scale: 75.0
L: 0.02

initial_conditions
T: 900.0
q: 0.0
A: 1.0
B: 0.0

stop_condition
enabled: true
species: A
threshold: 1e-4

solver_options
max_iterations: 42
rtol: 1e-6
atol: 1e-8
max_step: 1e-3
first_step: Some(1e-6)
vectorized: true
parallel: false
neighborhood_check: 1e-3
lsode2_method_family: auto
lsode2_symbolic_assembly: AtomView
lsode2_symbolic_execution: AOT
lsode2_aot_toolchain: c_tcc
lsode2_aot_profile: debug
lsode2_aot_output_dir: ./out
lsode2_linear_structure: banded
lsode2_banded_kl: 2
lsode2_banded_ku: 3
lsode2_linear_solver_policy: dense_lu
lsode2_native_execution: probe_before_bridge
lsode2_native_max_step_attempts: 111
lsode2_native_max_accepted_steps: 222
"#;

    let spec = parse_reactor_ivp_document_from_str(input)
        .expect("reactor IVP document should parse successfully");
    let serialized = spec.document_to_string();
    let roundtripped = parse_reactor_ivp_document_from_str(&serialized)
        .expect("serialized reactor IVP document should parse successfully");

    assert_eq!(
        roundtripped.problem_name.as_deref(),
        Some("RoundTripCondensed")
    );
    assert_eq!(
        roundtripped.problem_description.as_deref(),
        Some("Roundtrip document")
    );
    assert_eq!(roundtripped.physical.ro, 1800.0);
    assert_eq!(roundtripped.physical.Cp, 1200.0);
    assert_eq!(roundtripped.physical.Lambda, 0.35);
    assert_eq!(roundtripped.physical.m, 0.02);
    assert_eq!(roundtripped.physical.L, 0.02);
    assert_eq!(roundtripped.scaling.L, 0.02);
    assert_eq!(roundtripped.initial_conditions["T"], 900.0);
    assert_eq!(roundtripped.initial_conditions["q"], 0.0);
    assert_eq!(roundtripped.initial_conditions["A"], 1.0);
    assert_eq!(roundtripped.initial_conditions["B"], 0.0);
    assert_eq!(roundtripped.thermal_effects, vec![125000.0]);
    assert!(roundtripped.stop_condition.enabled);
    assert_eq!(roundtripped.stop_condition.species.as_deref(), Some("A"));
    assert_eq!(roundtripped.stop_condition.threshold, Some(1.0e-4));
    assert_eq!(
        roundtripped.solver_settings.solver.method,
        IvpMethodSpec::Lsode2
    );

    let lsode2 = roundtripped
        .solver_settings
        .solver_options
        .lsode2
        .as_ref()
        .expect("LSODE2-specific options should be preserved");
    assert!(matches!(
        lsode2.symbolic_execution,
        Some(Lsode2TaskExecutionSpec::Aot { .. })
    ));
    assert!(serialized.contains("stop_condition"));
    assert_eq!(serialized, roundtripped.document_to_string());
}
