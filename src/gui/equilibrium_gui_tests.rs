//! `egui_kittest` smoke coverage for the first equilibrium editor surface.

use super::equilibrium_gui::EquilibriumApp;
use super::equilibrium_gui_model::{
    CandidatePolicyDraft, EquilibriumGuiDocument, EquilibriumInventoryDraft,
    EquilibriumLookupDraft, EquilibriumProblemDraft, EquilibriumSolverDraft, GuiElementSearchMode,
    GuiPhSolveMode, GuiSolverBackend, GuiTraceSeedPolicyDraft,
};
use super::equilibrium_gui_request::{
    EquilibriumGuiSolveRequest, build_equilibrium_facade_request,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::Solvers;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverBackend, SolverPolicy,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::ChemEquilibrium::prelude::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition, PhaseComponentId, PhaseId,
    PhaseModel, PhaseSpec, PhysicalState, ResolvedThermochemistry, SubstanceSystemFactory,
    SubstanceSystemSpec, ThermoRepository,
};
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use crate::library_manager::with_library_manager;
use egui::accesskit::Role;
use egui_kittest::Harness;
use egui_kittest::kittest::Queryable;
use serde_json::json;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs;
use std::rc::Rc;
use std::sync::Arc;

fn wait_for_gui_worker(app: &mut EquilibriumApp) {
    for _ in 0..600 {
        app.poll_worker();
        if matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
                | super::equilibrium_gui_execution::EquilibriumGuiRunState::Failed { .. }
        ) {
            return;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    }
}

fn local_gui_library_snapshot() -> Vec<(String, Vec<u8>)> {
    let paths = with_library_manager(|manager| {
        vec![
            manager.substance_base_path().to_string(),
            manager.all_keys_substance_path().to_string(),
            manager.elements_path().to_string(),
        ]
    });
    paths
        .into_iter()
        .map(|path| {
            let bytes = fs::read(&path)
                .unwrap_or_else(|error| panic!("must read GUI fixture file '{path}': {error}"));
            (path, bytes)
        })
        .collect()
}

fn local_water_phase_app(
    temperature: super::equilibrium_gui_model::TemperatureDraft,
) -> EquilibriumApp {
    let mut app = EquilibriumApp::new();
    let mut gas_components = vec![super::equilibrium_gui_model::ComponentDraft {
        substance: "H2O".into(),
        initial_moles: "0.5".into(),
        source_library: Some("NASA_gas".into()),
    }];
    gas_components.push(super::equilibrium_gui_model::ComponentDraft {
        substance: "O2".into(),
        initial_moles: "0.25".into(),
        source_library: Some("NASA_gas".into()),
    });
    app.document.config.inventory =
        super::equilibrium_gui_model::EquilibriumInventoryDraft::ExplicitPhases {
            phases: vec![
                super::equilibrium_gui_model::PhaseDraft {
                    id: "gas".into(),
                    physical_state: super::equilibrium_gui_model::GuiPhysicalState::Gas,
                    model: super::equilibrium_gui_model::GuiPhaseModel::IdealGas,
                    components: gas_components,
                },
                super::equilibrium_gui_model::PhaseDraft {
                    id: "solid".into(),
                    physical_state: super::equilibrium_gui_model::GuiPhysicalState::Solid,
                    model: super::equilibrium_gui_model::GuiPhaseModel::PureCondensed,
                    components: vec![super::equilibrium_gui_model::ComponentDraft {
                        substance: "H2O(s)".into(),
                        initial_moles: "0".into(),
                        source_library: Some("NASA_cond".into()),
                    }],
                },
            ],
        };
    app.document.config.problem = super::equilibrium_gui_model::EquilibriumProblemDraft::FixedPt {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        temperature,
    };
    app.document.config.lookup = super::equilibrium_gui_model::EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into(), "NASA_cond".into()],
        permitted_libraries: vec!["NASA_gas".into(), "NASA_cond".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    app.document.config.diagnostics.collect_timing = true;
    app.document.config.phase_mode =
        super::equilibrium_gui_model::EquilibriumPhaseModeDraft::Bounded {
            // Match the engine's default temperature-scaled hysteresis at the
            // 250 K ice fixture. The GUI currently exposes explicit numeric
            // thresholds, so this keeps the story physically equivalent to
            // the canonical live-data test while still exercising the typed
            // GUI policy boundary.
            phase_epsilon: "1e-30".into(),
            dg_create: "-0.0020786".into(),
            dg_keep: "0.000020786".into(),
            max_phase_iterations: "20".into(),
            initial_phase_policy:
                super::equilibrium_gui_model::GuiInitialPhasePolicyDraft::AllDeclaredCandidates,
        };
    app
}

fn local_n2_fallback_app() -> EquilibriumApp {
    let mut app = EquilibriumApp::new();
    app.document.config.inventory = EquilibriumInventoryDraft::ExplicitPhases {
        phases: vec![super::equilibrium_gui_model::PhaseDraft {
            id: "gas".into(),
            physical_state: super::equilibrium_gui_model::GuiPhysicalState::Gas,
            model: super::equilibrium_gui_model::GuiPhaseModel::IdealGas,
            components: vec![
                super::equilibrium_gui_model::ComponentDraft {
                    substance: "N2".into(),
                    initial_moles: "1e-5".into(),
                    source_library: Some("NASA_gas".into()),
                },
                super::equilibrium_gui_model::ComponentDraft {
                    substance: "N".into(),
                    initial_moles: "1.0".into(),
                    source_library: Some("NASA_gas".into()),
                },
            ],
        }],
    };
    app.document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    app.document.config.problem = EquilibriumProblemDraft::FixedPt {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        temperature: super::equilibrium_gui_model::TemperatureDraft::Point {
            temperature_k: "5500".into(),
        },
    };
    app.document.config.solver.selection = EquilibriumSolverDraft::CustomCascade {
        backends: vec![GuiSolverBackend::RstNielsenLm, GuiSolverBackend::LegacyNr],
    };
    // Keep the first attempt deliberately short while leaving the legacy
    // fallback enough room to demonstrate its own step-control path.
    app.document.config.solver.overrides.max_iterations = "8".into();
    app.document.config.diagnostics.collect_timing = true;
    app
}

/// Builds the extensive enthalpy target used by the P,H N2/N GUI stories.
///
/// The target is derived through the same local repository and canonical P,T
/// solve as the worker request. This keeps the P,H stories independent from a
/// hard-coded number while preserving a stable real-data fixture.
fn local_h2o_ph_target() -> f64 {
    let data = ThermoData::try_new_fresh().expect("local catalog must load");
    let phase_id = PhaseId::new(Some("gas".into()));
    let phase = PhaseSpec::new(
        phase_id.clone(),
        vec!["H2".into(), "O2".into(), "H2O".into()],
        PhysicalState::Gas,
        PhaseModel::IdealGas,
    )
    .expect("H2O gas phase is valid");
    let spec = SubstanceSystemSpec::from_phases(vec![phase])
        .expect("H2O phase specification is valid")
        .with_lookup_policy(
            vec!["NASA_gas".into()],
            vec!["NASA_gas".into()],
            None,
            false,
        );
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        spec,
        Arc::clone(&data.repository),
    )
    .expect("H2O phase resolves");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("H2O layout builds");
    let composition = MultiphaseInitialComposition::from_sparse(
        &layout,
        vec![
            (PhaseComponentId::new(phase_id.clone(), "H2"), 0.1),
            (PhaseComponentId::new(phase_id.clone(), "O2"), 0.05),
            (PhaseComponentId::new(phase_id, "H2O"), 1.9),
        ],
    )
    .expect("H2O composition builds");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("H2O thermochemistry builds");
    let point = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(1_050.0, 101_325.0, 101_325.0)
                .expect("H2O P,T conditions are valid"),
            composition,
        )
        .with_phase_control_policy(PhaseControlPolicy::default()),
    )
    .expect("H2O reference P,T solve succeeds");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.component_moles(), 1_050.0)
        .expect("H2O reference enthalpy is finite")
}

fn configure_h2o_ph_app(
    target_enthalpy_j: f64,
    selection: EquilibriumSolverDraft,
) -> EquilibriumApp {
    let mut app = EquilibriumApp::new();
    app.document.config.solver.selection = selection;
    app.document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: target_enthalpy_j.to_string(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft {
            lower_k: "900".into(),
            upper_k: "1100".into(),
            seed_k: "1050".into(),
        },
    };
    if let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut app.document.config.inventory
    {
        phases[0].components = vec![
            super::equilibrium_gui_model::ComponentDraft {
                substance: "H2".into(),
                initial_moles: "0.1".into(),
                source_library: Some("NASA_gas".into()),
            },
            super::equilibrium_gui_model::ComponentDraft {
                substance: "O2".into(),
                initial_moles: "0.05".into(),
                source_library: Some("NASA_gas".into()),
            },
            super::equilibrium_gui_model::ComponentDraft {
                substance: "H2O".into(),
                initial_moles: "1.9".into(),
                source_library: Some("NASA_gas".into()),
            },
        ];
    }
    app.document.config.diagnostics.collect_timing = true;
    app
}

/// Derives a real water gas/ice P,H target at 250 K for the GUI transition
/// story. The reference solve uses the same bounded phase policy as the GUI
/// request, so the resulting target includes the phase-control inventory.
fn local_water_ph_target() -> f64 {
    let data = ThermoData::try_new_fresh().expect("local catalog must load");
    let gas_id = PhaseId::new(Some("gas".into()));
    let solid_id = PhaseId::new(Some("solid".into()));
    let gas = PhaseSpec::new(
        gas_id.clone(),
        vec!["H2O".into(), "O2".into()],
        PhysicalState::Gas,
        PhaseModel::IdealGas,
    )
    .expect("water gas phase is valid");
    let solid = PhaseSpec::new(
        solid_id.clone(),
        vec!["H2O(s)".into()],
        PhysicalState::Solid,
        PhaseModel::PureCondensed,
    )
    .expect("water solid phase is valid");
    let spec = SubstanceSystemSpec::from_phases(vec![gas, solid])
        .expect("water gas/solid specification is valid")
        .with_lookup_policy(
            vec!["NASA_gas".into(), "NASA_cond".into()],
            vec!["NASA_gas".into(), "NASA_cond".into()],
            None,
            false,
        );
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        spec,
        Arc::clone(&data.repository),
    )
    .expect("water gas/solid system resolves");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("water gas/solid layout builds");
    let composition = MultiphaseInitialComposition::from_sparse(
        &layout,
        vec![
            (PhaseComponentId::new(gas_id.clone(), "H2O"), 0.5),
            (PhaseComponentId::new(gas_id, "O2"), 0.25),
            (PhaseComponentId::new(solid_id, "H2O(s)"), 0.0),
        ],
    )
    .expect("water gas/solid composition builds");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("water gas/solid thermochemistry builds");
    let point = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(250.0, 101_325.0, 101_325.0)
                .expect("water P,T conditions are valid"),
            composition,
        )
        .with_phase_control_policy(
            PhaseControlPolicy::with_explicit_hysteresis(1e-30, -0.0020786, 0.000020786)
                .expect("water GUI hysteresis policy is valid"),
        ),
    )
    .expect("water gas/solid reference solve succeeds");
    thermochemistry
        .enthalpy_model()
        .evaluate_total(point.component_moles(), 250.0)
        .expect("water reference enthalpy is finite")
}

fn local_gas_continuation_app() -> EquilibriumApp {
    let mut app = EquilibriumApp::new();
    app.document.config.inventory = EquilibriumInventoryDraft::ExplicitPhases {
        phases: vec![super::equilibrium_gui_model::PhaseDraft {
            id: "gas".into(),
            physical_state: super::equilibrium_gui_model::GuiPhysicalState::Gas,
            model: super::equilibrium_gui_model::GuiPhaseModel::IdealGas,
            components: vec![
                super::equilibrium_gui_model::ComponentDraft {
                    substance: "H2O".into(),
                    initial_moles: "0.5".into(),
                    source_library: Some("NASA_gas".into()),
                },
                super::equilibrium_gui_model::ComponentDraft {
                    substance: "O2".into(),
                    initial_moles: "0.25".into(),
                    source_library: Some("NASA_gas".into()),
                },
            ],
        }],
    };
    app.document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    app.document.config.problem = EquilibriumProblemDraft::FixedPt {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        temperature: super::equilibrium_gui_model::TemperatureDraft::Range {
            start_k: "1000".into(),
            end_k: "1100".into(),
            point_count: "3".into(),
        },
    };
    app.document.config.solver.selection = EquilibriumSolverDraft::SingleBackend {
        backend: GuiSolverBackend::LegacyNr,
    };
    app.document.config.diagnostics.collect_timing = true;
    app.document.config.postprocessing.plot_target =
        super::equilibrium_gui_model::GuiPlotTarget::Embedded;
    app
}

fn assert_result_diagnostics_are_rendered(app: EquilibriumApp, first_point_label: &str) {
    let app = Rc::new(RefCell::new(app));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Results")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, first_point_label)
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solve diagnostics")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Conservation");
    harness.get_by_label("Residuals");
    harness.get_by_label("Fallback attempts");
    harness.get_by_label("K_eq validation");
    harness.get_by_label("Lookup provenance");
}

fn assert_phase_lifecycle_trace_is_rendered(app: EquilibriumApp, first_point_label: &str) {
    let app = Rc::new(RefCell::new(app));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Results")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, first_point_label)
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase lifecycle")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Decision");
    harness.get_by_label("Evidence");
}

fn assert_ph_result_diagnostics_are_rendered(app: EquilibriumApp) {
    let app = Rc::new(RefCell::new(app));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Results")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "P,H solve diagnostics")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Inner backend attempts");
    harness.get_by_label("Phase-control transitions");
    harness.get_by_label("Outer solve path");
}

#[test]
fn editor_exposes_canonical_request_and_validation_controls() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_role_and_label(Role::Button, "Validate document");
    harness.get_by_role_and_label(Role::Button, "Prepare canonical request");
    harness
        .get_by_role_and_label(Role::Button, "Libraries")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Library lookup")
        .click_accesskit();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Diagnostics")
        .click_accesskit();
    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Equilibrium-constant validation");
}

#[test]
fn lookup_and_diagnostics_controls_render_typed_policies() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Libraries")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Library lookup")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Engine default");
    harness.get_by_label("Local catalog: not loaded");
    harness
        .get_by_role_and_label(Role::Button, "Declared component lookup")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Lookup instruction");
    harness.get_by_label("engine policy");
    harness.get_by_label("Explicit policy").click_accesskit();
    harness.run();
    harness.get_by_label("Priority libraries (ordered)");
    harness.get_by_label("Permitted libraries (closed candidate set)");
    harness.get_by_label("Allow online NIST fallback");
    assert!(matches!(
        &app.borrow().document.config.lookup,
        EquilibriumLookupDraft::Explicit { .. }
    ));

    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Diagnostics")
        .click_accesskit();
    harness.run();
    harness.get_by_role_and_label(Role::CheckBox, "Collect timing");
    assert!(harness.query_by_label("Retain backend attempts").is_none());
    assert!(
        harness
            .query_by_label("Retain conservation report")
            .is_none()
    );
    assert!(harness.query_by_label("Retain phase transitions").is_none());
    assert!(harness.query_by_label("Phase lifecycle trace").is_none());
    assert!(
        harness
            .query_by_label("Retained lifecycle events override")
            .is_none()
    );
    assert!(harness.query_by_label("Range lifecycle trace").is_none());
    assert!(!app.borrow().document.config.diagnostics.collect_timing);
    assert_eq!(
        app.borrow()
            .document
            .config
            .diagnostics
            .phase_lifecycle_trace,
        super::equilibrium_gui_model::GuiPhaseLifecycleTrace::Off
    );

    app.borrow_mut().document.config.phase_mode =
        super::equilibrium_gui_model::EquilibriumPhaseModeDraft::Bounded {
            phase_epsilon: "1e-12".into(),
            dg_create: "-1e-6".into(),
            dg_keep: "1e-8".into(),
            max_phase_iterations: "20".into(),
            initial_phase_policy:
                super::equilibrium_gui_model::GuiInitialPhasePolicyDraft::FromInitialMoles,
        };
    if let EquilibriumProblemDraft::FixedPt { temperature, .. } =
        &mut app.borrow_mut().document.config.problem
    {
        *temperature = super::equilibrium_gui_model::TemperatureDraft::Range {
            start_k: "300".into(),
            end_k: "1000".into(),
            point_count: "5".into(),
        };
    }
    harness.run();
    harness.get_by_label("Phase lifecycle trace");
    harness.get_by_label("Retained lifecycle events override");
    harness.get_by_label("Range lifecycle trace");
}

#[test]
fn element_mode_exposes_candidate_preview_controls_in_egui() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut().document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["C".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft::default(),
        assignments: vec![Default::default()],
    };
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Setup")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Search by elements");
    harness.get_by_label("Candidate preview is required before phase assignment");
    harness.get_by_role_and_label(Role::Button, "Preview candidates");
    harness.get_by_label("Exact set");
    harness.get_by_label("Max candidates");
}

#[test]
fn phase_editor_keeps_semantic_row_labels_when_phases_are_reordered() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut().document.config.inventory = EquilibriumInventoryDraft::ExplicitPhases {
        phases: vec![
            super::equilibrium_gui_model::PhaseDraft {
                id: "gas".into(),
                physical_state: super::equilibrium_gui_model::GuiPhysicalState::Gas,
                model: super::equilibrium_gui_model::GuiPhaseModel::IdealGas,
                components: vec![Default::default()],
            },
            super::equilibrium_gui_model::PhaseDraft {
                id: "liquid".into(),
                physical_state: super::equilibrium_gui_model::GuiPhysicalState::Liquid,
                model: super::equilibrium_gui_model::GuiPhaseModel::PureCondensed,
                components: vec![Default::default()],
            },
        ],
    };
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Setup")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Phase 'gas'");
    harness.get_by_label("Phase 'liquid'");

    // The rows are scoped by semantic phase IDs, not by their vector index.
    // Reordering the document must therefore preserve the same two editor
    // identities and only change their visual order.
    if let EquilibriumInventoryDraft::ExplicitPhases { phases } =
        &mut app.borrow_mut().document.config.inventory
    {
        phases.swap(0, 1);
    }
    harness.run();
    harness.get_by_label("Phase 'liquid'");
    harness.get_by_label("Phase 'gas'");
}

#[test]
fn fixed_ph_controls_are_visible_in_egui() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "0".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("P,H = const");
    harness.get_by_label("Target total enthalpy [J]");
    harness.get_by_label("Lower bound [K]");
    harness.get_by_label("Initial seed [K]");
}

#[test]
fn solver_exposes_ph_route_only_for_a_ph_problem() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "0".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness.get_by_label("P,H route");
    assert_eq!(
        app.borrow().document.config.solver.ph_solve_mode,
        GuiPhSolveMode::Auto
    );

    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::default();
    harness.run();
    assert!(harness.query_by_label("P,H route").is_none());
}

#[test]
fn problem_tab_shows_only_controls_for_the_selected_thermodynamic_route() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("P,T = const");
    harness.get_by_label("Point");
    assert!(
        harness
            .query_by_label("Target total enthalpy [J]")
            .is_none()
    );

    app.borrow_mut().document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    harness.run();
    harness.get_by_label("P,H = const");
    harness.get_by_label("Target total enthalpy [J]");
    harness.get_by_label("Lower bound [K]");
    assert!(harness.query_by_label("Point").is_none());
}

#[test]
fn fixed_ph_document_changes_mark_the_prepared_request_stale() {
    let mut app = EquilibriumApp::new();
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "1000".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    app.prepare_request().expect("valid P,H request prepares");
    assert!(app.prepared_request_is_current());

    // The editor owns the fingerprint boundary. A P,H target edit must drop
    // the prepared request just like a P,T or inventory edit, otherwise the
    // worker could solve with an old energy target.
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "2000".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    // Direct model mutation is intentionally observable as stale even before
    // the next editor frame clears the derived request slot. The Run control
    // is gated by this predicate, so an old request cannot execute.
    assert!(app.prepared_request().is_some());
    assert!(!app.prepared_request_is_current());
}

#[test]
fn solver_panel_exposes_trace_seed_override_without_changing_defaults() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    app.borrow_mut()
        .document
        .config
        .solver
        .overrides
        .trace_seed_policy = Some(GuiTraceSeedPolicyDraft::RelativeToLargestInitialMole {
        fraction: "1e-12".into(),
        minimum_floor: "1e-30".into(),
    });
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Override trace-species seed policy");
    harness.get_by_label("Trace seed strategy");
    harness.get_by_label("Relative to largest mole");
    harness.get_by_label("Trace fraction");
    harness.get_by_label("Minimum trace floor");
}

#[test]
fn solver_editor_exposes_optional_cascade_budget_without_changing_defaults() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    assert!(
        app.borrow()
            .document
            .config
            .solver
            .overrides
            .cascade_budget
            .is_none()
    );
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness
        .get_by_label("Override cascade budget")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Cascade max attempts");
    harness.get_by_label("Cascade iterations per attempt");
    harness.get_by_label("Cascade total iterations");
    assert!(
        app.borrow()
            .document
            .config
            .solver
            .overrides
            .cascade_budget
            .is_some()
    );
}

#[test]
fn custom_solver_cascade_is_ordered_and_rejects_duplicates() {
    let mut document = EquilibriumGuiDocument::new();
    document.config.solver.selection = EquilibriumSolverDraft::CustomCascade {
        backends: vec![GuiSolverBackend::LegacyNr, GuiSolverBackend::RstLm],
    };
    let validated = document
        .validate_for_run()
        .expect("a unique custom cascade is valid");
    let request =
        build_equilibrium_facade_request(validated, None).expect("cascade facade request builds");
    assert!(matches!(request, EquilibriumGuiSolveRequest::Facade(_)));

    document.config.solver.selection = EquilibriumSolverDraft::CustomCascade {
        backends: vec![GuiSolverBackend::LegacyNr, GuiSolverBackend::LegacyNr],
    };
    let report = document
        .validate()
        .expect_err("duplicate cascade entries must be rejected");
    assert!(
        report
            .issues
            .iter()
            .any(|issue| issue.field == "solver.selection.backends[1]")
    );
}

#[test]
fn solver_editor_exposes_custom_cascade_controls() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Numerics")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Solver")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Use custom cascade").click_accesskit();
    harness.run();
    harness.get_by_label("Ordered fallback sequence");
    harness.get_by_label("Add backend");
}

#[test]
fn validation_errors_are_visible_after_editing_a_numeric_field() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    if let super::equilibrium_gui_model::EquilibriumProblemDraft::FixedPt { pressure_pa, .. } =
        &mut app.borrow_mut().document.config.problem
    {
        *pressure_pa = "0".into();
    }
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Validate document")
        .click_accesskit();
    harness.run();
    harness.get_by_label("problem.pressure_pa: value must be greater than zero");
}

#[test]
fn every_concrete_backend_can_be_prepared_through_the_gui_boundary() {
    for backend in GuiSolverBackend::ALL {
        let mut app = EquilibriumApp::new();
        app.document.config.solver.selection =
            super::equilibrium_gui_model::EquilibriumSolverDraft::SingleBackend { backend };
        app.prepare_request()
            .unwrap_or_else(|error| panic!("{} must prepare: {error}", backend.label()));
    }
}

#[test]
fn prepare_button_reports_a_valid_default_document() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    assert!(matches!(
        app.borrow().document.config.lookup,
        EquilibriumLookupDraft::Default
    ));
    assert!(matches!(
        &app.borrow().document.config.solver.selection,
        EquilibriumSolverDraft::ProductionDefault
    ));
    assert!(!app.borrow().document.config.diagnostics.collect_timing);
    assert!(matches!(
        &app.borrow()
            .document
            .config
            .diagnostics
            .phase_lifecycle_trace,
        super::equilibrium_gui_model::GuiPhaseLifecycleTrace::Off
    ));
    assert!(matches!(
        &app.borrow().document.config.postprocessing.plot_target,
        super::equilibrium_gui_model::GuiPlotTarget::None
    ));
    assert!(matches!(
        &app.borrow().document.config.postprocessing.resampling,
        super::equilibrium_gui_model::GuiResamplingDraft::None
    ));
    assert!(
        !app.borrow()
            .document
            .config
            .solver
            .overrides
            .scaling_enabled,
        "GUI default must preserve the engine's stable scaling default"
    );
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Prepare canonical request")
        .click_accesskit();
    harness.run();
    assert!(app.borrow().prepared_request().is_some());
    assert!(app.borrow().prepared_request_is_current());
    assert!(matches!(
        app.borrow().prepared_request(),
        Some(EquilibriumGuiSolveRequest::Facade(_))
    ));
    harness.get_by_role_and_label(Role::Button, "Run prepared request");
}

#[test]
fn phase_policy_exposes_bounded_hysteresis_controls() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Phase policy")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Fixed declared phases");
    assert!(harness.query_by_label("Phase epsilon").is_none());
    harness
        .get_by_label("Bounded phase control")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Phase epsilon");
    harness.get_by_label("Creation driving force");
    harness.get_by_label("Keep driving force");
    harness.get_by_label("Maximum phase iterations");
    harness.get_by_label("Initial phase set");
    harness
        .get_by_label("Fixed declared phases")
        .click_accesskit();
    harness.run();
    assert!(harness.query_by_label("Phase epsilon").is_none());
    assert!(harness.query_by_label("Initial phase set").is_none());
}

#[test]
fn edited_document_cannot_publish_an_old_worker_result() {
    let mut app = EquilibriumApp::new();
    app.prepare_request().expect("default request prepares");
    let ticket = app.begin_prepared_run().expect("run ticket exists");
    if let super::equilibrium_gui_model::EquilibriumProblemDraft::FixedPt {
        temperature: super::equilibrium_gui_model::TemperatureDraft::Point { temperature_k },
        ..
    } = &mut app.document.config.problem
    {
        *temperature_k = "1100".into();
    } else {
        panic!("default document must use a point temperature");
    }

    let publication = app.publish_failure(ticket, "late worker failure");
    assert_eq!(
        publication,
        super::equilibrium_gui_execution::EquilibriumGuiPublication::Stale
    );
    assert!(app.result_snapshot().is_none());
}

#[test]
fn failed_current_worker_keeps_previous_result_slot_empty_and_exposes_error() {
    let mut app = EquilibriumApp::new();
    app.prepare_request().expect("default request prepares");
    let ticket = app.begin_prepared_run().expect("run ticket exists");
    assert_eq!(
        app.publish_failure(ticket, "solver worker failed"),
        super::equilibrium_gui_execution::EquilibriumGuiPublication::Accepted
    );
    assert!(app.result_snapshot().is_none());
    assert_eq!(app.last_error(), Some("solver worker failed"));
}

#[test]
fn app_lifecycle_gate_can_cancel_before_a_worker_is_started() {
    let mut app = EquilibriumApp::new();
    app.prepare_request().expect("default request prepares");
    let ticket = app.begin_prepared_run().expect("run ticket exists");
    assert_eq!(
        app.run_state(),
        super::equilibrium_gui_execution::EquilibriumGuiRunState::Solving {
            run_id: ticket.run_id()
        }
    );
    assert!(app.cancel_run());
    assert!(app.result_snapshot().is_none());
}

#[test]
fn plot_visibility_is_display_state_and_does_not_stale_prepared_request() {
    let mut app = EquilibriumApp::new();
    app.prepare_request().expect("default request prepares");
    assert!(app.prepared_request_is_current());

    app.set_plot_series_visible("gas::H2O", false);

    assert!(!app.plot_series_visible("gas::H2O"));
    assert!(app.prepared_request_is_current());
    app.set_plot_series_visible("gas::H2O", true);
    assert!(app.plot_series_visible("gas::H2O"));
}

#[test]
fn range_editor_exposes_grid_controls_and_accepts_descending_input() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_label("Range").click_accesskit();
    harness.run();
    harness.get_by_label("Start [K]");
    harness.get_by_label("End [K]");
    harness.get_by_label("Solved points");

    let mut app = app.borrow_mut();
    if let super::equilibrium_gui_model::EquilibriumProblemDraft::FixedPt {
        temperature:
            super::equilibrium_gui_model::TemperatureDraft::Range {
                start_k,
                end_k,
                point_count,
            },
        ..
    } = &mut app.document.config.problem
    {
        *start_k = "1500".into();
        *end_k = "300".into();
        *point_count = "4".into();
    } else {
        panic!("range control must create a range draft");
    }
    assert!(app.prepare_request().is_ok());
}

#[test]
fn output_shows_resampling_only_for_a_pt_temperature_range() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Output")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Postprocessing and plots")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Result table");
    assert!(harness.query_by_label("PCHIP display resampling").is_none());

    if let EquilibriumProblemDraft::FixedPt { temperature, .. } =
        &mut app.borrow_mut().document.config.problem
    {
        *temperature = super::equilibrium_gui_model::TemperatureDraft::Range {
            start_k: "300".into(),
            end_k: "1000".into(),
            point_count: "5".into(),
        };
    } else {
        panic!("the default GUI problem must be P,T");
    }
    harness.run();
    harness.get_by_label("PCHIP display resampling");
}

#[test]
fn postprocessing_changes_do_not_invalidate_a_prepared_solver_request() {
    let mut app = EquilibriumApp::new();
    app.prepare_request().expect("default request prepares");
    app.document.config.postprocessing.result_basis =
        super::equilibrium_gui_model::GuiResultBasis::PhaseTotals;
    app.document.config.postprocessing.resampling =
        super::equilibrium_gui_model::GuiResamplingDraft::Pchip {
            output_points: "200".into(),
            interpolation_space: super::equilibrium_gui_model::GuiInterpolationSpace::Log,
            clamp: true,
        };
    app.document.config.postprocessing.y_scale = super::equilibrium_gui_model::GuiPlotScale::Log10;
    app.document.config.postprocessing.table_density =
        super::equilibrium_gui_model::GuiResultTableDensity::Compact;
    assert!(app.prepared_request().is_some());
    assert!(app.prepared_request_is_current());
}

#[test]
#[ignore = "requires the local thermochemical catalog and runs a real solver worker"]
fn offline_local_h2o_point_story_publishes_provenance() {
    let mut app = EquilibriumApp::new();
    // The default GUI policy is intentionally repository-owned and may not
    // expose every local catalog alias. Pin this real fixture explicitly so
    // the story tests the solve/lifecycle contract, not catalog ordering.
    app.document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    if let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut app.document.config.inventory
    {
        phases[0]
            .components
            .push(super::equilibrium_gui_model::ComponentDraft {
                substance: "O2".into(),
                initial_moles: "0.25".into(),
                source_library: Some("NASA_gas".into()),
            });
    }
    app.prepare_request()
        .expect("the explicit local H2O document prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert_eq!(
        app.run_state(),
        super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { run_id: 1 },
        "default local H2O solve failed: {:?}",
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("successful local solve publishes a snapshot");
    assert_eq!(snapshot.points().len(), 1);
    assert!(!snapshot.component_labels().is_empty());
    assert!(
        !snapshot.points()[0]
            .source()
            .build_report()
            .components()
            .is_empty()
    );
    assert!(snapshot.points()[0].source().solve_report().attempt_count() >= 1);
    assert!(
        snapshot.points()[0]
            .source()
            .summary_rows()
            .iter()
            .any(|row| { matches!(row.section, "backend" | "validation" | "acceptance") })
    );
    assert_result_diagnostics_are_rendered(app, "Point 1: 1000.000000 K");
}

#[test]
#[ignore = "requires the local NASA gas catalog and runs a real P,H GUI worker"]
fn offline_local_h2o_ph_story_publishes_energy_contract() {
    let data = ThermoData::try_new_fresh().expect("local catalog must load");
    let phase_id = PhaseId::new(Some("gas".into()));
    let phase = PhaseSpec::new(
        phase_id.clone(),
        vec!["H2".into(), "O2".into(), "H2O".into()],
        PhysicalState::Gas,
        PhaseModel::IdealGas,
    )
    .expect("H2O gas phase is valid");
    let spec = SubstanceSystemSpec::from_phases(vec![phase])
        .expect("phase specification is valid")
        .with_lookup_policy(
            vec!["NASA_gas".into()],
            vec!["NASA_gas".into()],
            None,
            false,
        );
    let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
        spec,
        Arc::clone(&data.repository),
    )
    .expect("live H2O phase resolves");
    let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
        .expect("resolved layout builds");
    let composition = MultiphaseInitialComposition::from_sparse(
        &layout,
        vec![
            (PhaseComponentId::new(phase_id.clone(), "H2"), 0.1),
            (PhaseComponentId::new(phase_id.clone(), "O2"), 0.05),
            (PhaseComponentId::new(phase_id, "H2O"), 1.9),
        ],
    )
    .expect("initial composition builds");
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
        .expect("live H2O thermochemistry bundle builds");
    let pt = solve_resolved_pt(
        ResolvedPhaseEquilibriumRequest::new(
            &resolved,
            EquilibriumConditions::new(2_500.0, 101_325.0, 101_325.0)
                .expect("live P,T conditions are valid"),
            composition.clone(),
        )
        .with_solve_options(
            crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::EquilibriumSolveOptions::new()
                .with_solver_policy(
                    SolverPolicy::Single(SolverBackend::Legacy(Solvers::NR)),
                )
                .expect("single legacy reference policy is valid"),
        ),
    )
    .expect("live reference P,T solve must be accepted");
    let target = thermochemistry
        .enthalpy_model()
        .evaluate_total(pt.component_moles(), 2_500.0)
        .expect("live H2O enthalpy evaluates");

    let mut app = EquilibriumApp::new();
    // Use one deterministic legacy backend here so this ignored GUI story
    // validates the route-dependent monolithic snapshot itself rather than
    // conflating it with the broader RST backend cascade characterization.
    app.document.config.solver.selection = EquilibriumSolverDraft::SingleBackend {
        backend: GuiSolverBackend::LegacyNr,
    };
    app.document.config.lookup = EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: target.to_string(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft {
            lower_k: "1900".into(),
            upper_k: "2900".into(),
            seed_k: "1900".into(),
        },
    };
    if let EquilibriumInventoryDraft::ExplicitPhases { phases } = &mut app.document.config.inventory
    {
        phases[0].components = vec![
            super::equilibrium_gui_model::ComponentDraft {
                substance: "H2".into(),
                initial_moles: "0.1".into(),
                source_library: Some("NASA_gas".into()),
            },
            super::equilibrium_gui_model::ComponentDraft {
                substance: "O2".into(),
                initial_moles: "0.05".into(),
                source_library: Some("NASA_gas".into()),
            },
            super::equilibrium_gui_model::ComponentDraft {
                substance: "H2O".into(),
                initial_moles: "1.9".into(),
                source_library: Some("NASA_gas".into()),
            },
        ];
    }
    app.prepare_request().expect("P,H GUI request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);
    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "live P,H GUI solve failed: {:?}",
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("P,H worker publishes an immutable snapshot");
    let energy = snapshot
        .enthalpy()
        .expect("P,H snapshot retains the energy contract");
    assert!((energy.calculated_enthalpy_j() - target).abs() < 1e-4);
    assert_eq!(snapshot.points().len(), 1);
    let ph_diagnostics = snapshot
        .ph_diagnostics()
        .expect("P,H snapshot retains outer-solver diagnostics");
    // The GUI uses the canonical coupled route. It must expose monolithic
    // evidence instead of manufacturing scalar temperature trials. Keep the
    // nested assertion as a compatibility guard for an explicitly selected
    // reference route.
    if let Some(monolithic) = ph_diagnostics.monolithic() {
        assert!(!monolithic.backend_attempts().is_empty());
        assert!(!monolithic.accepted_backend().is_empty());
        // Legacy fallback metrics may intentionally omit per-evaluation
        // counters; backend-attempt and acceptance evidence remain required.
    } else {
        assert!(ph_diagnostics.trial_count() >= 1);
    }
    assert!(ph_diagnostics.solved_temperature_k().is_finite());
    assert!(!ph_diagnostics.solve_path().is_empty());
    assert!(ph_diagnostics.inner_backend_attempts() >= 1);

    let app = Rc::new(RefCell::new(app));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "Results")
        .click_accesskit();
    harness.run();
    harness
        .get_by_role_and_label(Role::Button, "P,H solve diagnostics")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Solved temperature [K]");
    harness.get_by_label("Outer solve path");
    if ph_diagnostics.monolithic().is_some() {
        harness
            .get_by_role_and_label(Role::Button, "P,H monolithic evidence")
            .click_accesskit();
        harness.run();
        harness.get_by_label("Accepted backend");
        harness.get_by_role_and_label(Role::Button, "Backend attempts");
    } else {
        harness
            .get_by_role_and_label(Role::Button, "P,H temperature trials")
            .click_accesskit();
        harness.run();
        harness.get_by_label("Step");
        harness.get_by_label("Scaled error");
    }
}

#[test]
#[ignore = "requires the local NASA gas catalog and exercises the P,H inner fallback"]
fn offline_local_h2o_ph_inner_fallback_story_renders_backend_attempts() {
    let target = local_h2o_ph_target();
    let mut app = configure_h2o_ph_app(
        target,
        EquilibriumSolverDraft::CustomCascade {
            // This narrow bracket forces the symbolic route to fall back to
            // the safeguarded numeric route; Nielsen is retained as the first
            // inner candidate and Legacy NR receives the recovery budget.
            backends: vec![GuiSolverBackend::RstNielsenLm, GuiSolverBackend::LegacyNr],
        },
    );
    app.document.config.solver.ph_solve_mode = GuiPhSolveMode::NestedTemperature;
    app.document.config.solver.overrides.max_iterations = "1000".into();
    app.prepare_request().expect("H2O P,H request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "P,H inner fallback did not complete: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("P,H fallback publishes an immutable snapshot");
    let diagnostics = snapshot
        .ph_diagnostics()
        .expect("P,H fallback retains route diagnostics");
    let attempt_count = diagnostics
        .monolithic()
        .map(|evidence| evidence.backend_attempts().len())
        .unwrap_or_else(|| diagnostics.inner_backend_attempts());
    assert!(
        attempt_count >= 2,
        "P,H inner fallback must retain both backend attempts: {diagnostics:?}"
    );
    assert!(diagnostics.fallback_reason().is_none());
}

#[test]
#[ignore = "requires the local NASA gas catalog and exercises P,H all-backends-failed"]
fn offline_local_h2o_ph_all_backends_failed_is_transactional() {
    let target = local_h2o_ph_target();
    let mut app = configure_h2o_ph_app(
        target,
        EquilibriumSolverDraft::SingleBackend {
            backend: GuiSolverBackend::RstNielsenLm,
        },
    );
    app.document.config.solver.ph_solve_mode = GuiPhSolveMode::NestedTemperature;
    // One iteration is deliberately insufficient for this real problem. The
    // resulting failure must remain an engine failure, not a GUI validation
    // error or a fabricated partial P,H result.
    app.document.config.solver.overrides.max_iterations = "1".into();
    app.prepare_request()
        .expect("the deliberately under-budget P,H request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Failed { .. }
        ),
        "under-budget P,H solve must fail visibly: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    assert!(
        app.last_error()
            .is_some_and(|error| error.contains("all equilibrium solver backends failed")),
        "failure must preserve the typed all-backends-failed diagnostic: {:?}",
        app.last_error()
    );
    assert!(app.result_snapshot().is_none());
}

#[test]
fn fixed_ph_cancellation_discards_late_worker_completion() {
    let mut app = EquilibriumApp::new();
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "1.0".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    app.prepare_request()
        .expect("P,H document must prepare before cancellation");
    let ticket = app
        .begin_prepared_run()
        .expect("P,H run ticket must be created");
    assert!(app.cancel_run());
    assert!(matches!(
        app.run_state(),
        super::equilibrium_gui_execution::EquilibriumGuiRunState::Cancelling { .. }
    ));
    assert!(app.result_snapshot().is_none());
    assert_eq!(
        app.publish_failure(ticket, "late P,H worker failure"),
        super::equilibrium_gui_execution::EquilibriumGuiPublication::Stale
    );
}

#[test]
#[ignore = "requires the local NASA gas/condensed catalogs and exercises P,H phase publication"]
fn offline_local_water_ph_story_publishes_phase_transition() {
    let before = local_gui_library_snapshot();
    let target = local_water_ph_target();
    let mut app = local_water_phase_app(super::equilibrium_gui_model::TemperatureDraft::Point {
        temperature_k: "250".into(),
    });
    // The mixed gas/condensed fixture intentionally exercises the legacy
    // numeric route: the symbolic monolithic P,H payload requires one native
    // coefficient interval for every component, while the GUI must still
    // publish a valid phase-transition result for this supported problem.
    app.document.config.solver.selection = EquilibriumSolverDraft::ProductionDefault;
    app.document.config.solver.overrides.max_iterations = "1000".into();
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: target.to_string(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft {
            lower_k: "240".into(),
            upper_k: "260".into(),
            seed_k: "250".into(),
        },
    };
    app.prepare_request()
        .expect("water gas/ice P,H request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "water gas/ice P,H solve failed: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("water gas/ice P,H result publishes");
    assert_eq!(snapshot.phase_labels(), &["gas", "solid"]);
    assert!(snapshot.points()[0].phase_totals()[1] > 0.49);
    let diagnostics = snapshot
        .ph_diagnostics()
        .expect("water P,H result retains route diagnostics");
    assert!(diagnostics.phase_control_transitions() >= 1);
    assert_eq!(before, local_gui_library_snapshot());
    assert_ph_result_diagnostics_are_rendered(app);
}

#[test]
#[ignore = "requires the local NASA gas/condensed catalogs and runs an unreachable P,H worker"]
fn offline_local_unreachable_ph_target_is_transactional() {
    let before = local_gui_library_snapshot();
    let mut app = local_water_phase_app(super::equilibrium_gui_model::TemperatureDraft::Point {
        temperature_k: "350".into(),
    });
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "101325".into(),
        reference_pressure_pa: "101325".into(),
        // Finite but far outside the resolved water enthalpy range. The
        // request is structurally valid; the worker must reject it without
        // publishing a fabricated endpoint or partial result.
        target_enthalpy_j: "1e99".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    app.prepare_request()
        .expect("an unreachable finite P,H target remains structurally valid");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Failed { .. }
        ),
        "unreachable P,H target must fail visibly: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    assert!(app.last_error().is_some());
    assert!(
        app.result_snapshot().is_none(),
        "an unreachable P,H target must not publish a partial snapshot"
    );
    assert_eq!(before, local_gui_library_snapshot());
}

#[test]
#[ignore = "requires the local NASA gas catalog and exercises a real fallback worker"]
fn offline_local_fallback_story_renders_backend_attempts() {
    let mut app = local_n2_fallback_app();
    app.prepare_request()
        .expect("the real N2/N fallback request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "fallback story did not complete: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("fallback solve publishes an immutable snapshot");
    let report = snapshot.points()[0].source().solve_report();
    // The facade preserves the declared cascade and always publishes the
    // accepted backend evidence. A fallback is a numerical outcome, not a UI
    // contract: current local coefficients allow Nielsen LM to accept before
    // Legacy NR is needed.
    assert!(
        report.attempt_count() >= 1,
        "missing backend evidence: {}",
        report.summary()
    );
    assert_result_diagnostics_are_rendered(app, "Point 1: 5500.000000 K");
}

#[test]
#[ignore = "requires the local NASA gas catalog and runs a real continuation worker"]
fn offline_local_gas_temperature_range_story_reuses_without_phase_transition() {
    let before = local_gui_library_snapshot();
    let mut app = local_gas_continuation_app();
    app.prepare_request()
        .expect("the real gas continuation request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "gas continuation did not complete: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("gas continuation publishes a range snapshot");
    assert_eq!(snapshot.points().len(), 3);
    assert_eq!(snapshot.phase_labels(), &["gas"]);
    let report = snapshot
        .range_report()
        .expect("range result retains continuation evidence");
    assert_eq!(report.formulation_builds(), 1);
    assert_eq!(report.formulation_reuses(), 2);
    assert_eq!(report.phase_control_transitions(), 0);
    assert!(report.point_timing().total() > std::time::Duration::ZERO);
    assert_eq!(before, local_gui_library_snapshot());
    app.open_embedded_plot()
        .expect("an accepted range snapshot must build an embedded plot without another solve");
    assert_result_diagnostics_are_rendered(app, "Point 1: 1000.000000 K");
}

#[test]
#[ignore = "requires the local NASA gas/condensed catalogs and runs a real multiphase worker"]
fn offline_local_water_ice_story_publishes_phase_totals_and_status() {
    let mut app = local_water_phase_app(super::equilibrium_gui_model::TemperatureDraft::Point {
        temperature_k: "250".into(),
    });
    app.document.config.diagnostics.phase_lifecycle_trace =
        super::equilibrium_gui_model::GuiPhaseLifecycleTrace::PhaseLifecycle;
    app.prepare_request()
        .expect("the real water/ice GUI request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);
    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "water/ice GUI solve did not complete: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app
        .result_snapshot()
        .expect("water/ice result is published");
    assert_eq!(snapshot.points().len(), 1);
    assert_eq!(snapshot.phase_labels(), &["gas", "solid"]);
    assert!(snapshot.points()[0].phase_totals()[1] > 0.49);
    assert!(
        snapshot.points()[0]
            .phase_statuses()
            .iter()
            .any(|status| status == "Active" || status == "Appeared")
    );
    assert!(snapshot.points()[0].source().solve_report().attempt_count() >= 1);
    let trace = snapshot.points()[0]
        .lifecycle_trace()
        .expect("enabled GUI trace is retained on the accepted point");
    assert!(
        trace
            .events()
            .iter()
            .any(|event| event.title().contains("TPD") || event.title().contains("transition")),
        "phase fixture must retain a physical stability or transition decision: {trace:?}"
    );
    assert_phase_lifecycle_trace_is_rendered(app, "Point 1: 250.000000 K");
}

#[test]
#[ignore = "requires the local NASA gas/condensed catalogs and runs a real range worker"]
fn offline_local_water_temperature_range_story_keeps_phase_layout() {
    let mut app = local_water_phase_app(super::equilibrium_gui_model::TemperatureDraft::Range {
        start_k: "250".into(),
        // The bundled NASA_cond H2O(s) record ends at 273.15 K. Keep this
        // GUI continuation story inside that real coefficient interval; a
        // phase-transition range belongs to the dedicated engine fixtures
        // with per-phase temperature policies.
        end_k: "270".into(),
        point_count: "3".into(),
    });
    app.prepare_request()
        .expect("the real water range GUI request prepares");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);
    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "water range GUI solve did not complete: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    let snapshot = app.result_snapshot().expect("range result is published");
    assert_eq!(snapshot.points().len(), 3);
    assert_eq!(snapshot.phase_labels(), &["gas", "solid"]);
    assert!(
        snapshot
            .points()
            .iter()
            .all(|point| point.phase_totals().len() == 2)
    );
    let report = snapshot
        .range_report()
        .expect("range snapshot retains continuation evidence");
    assert_eq!(report.point_count(), 3);
    assert_eq!(report.formulation_builds(), 1);
    assert_eq!(report.formulation_reuses(), 2);
    assert!(
        report.phase_control_transitions() > 0,
        "real water range must retain at least one phase transition"
    );
    assert!(report.point_timing().total() > std::time::Duration::ZERO);
    assert_result_diagnostics_are_rendered(app, "Point 1: 250.000000 K");
}

#[test]
#[ignore = "requires the local NASA gas/condensed catalogs and runs a failing range worker"]
fn offline_local_water_temperature_range_failure_is_transactional() {
    let before = local_gui_library_snapshot();
    let mut app = local_water_phase_app(super::equilibrium_gui_model::TemperatureDraft::Range {
        start_k: "250".into(),
        // H2O(s) in the bundled NASA_cond fixture ends at 273.15 K. This
        // deliberately crosses that boundary so the worker must report a
        // typed range failure instead of publishing early points.
        end_k: "300".into(),
        point_count: "3".into(),
    });
    app.prepare_request()
        .expect("the failing range remains a structurally valid request");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);

    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Failed { .. }
        ),
        "out-of-coverage range must fail visibly: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    assert!(app.last_error().is_some());
    assert!(
        app.result_snapshot().is_none(),
        "a failed range must not publish a partial result"
    );
    assert_eq!(before, local_gui_library_snapshot());
}

#[test]
fn document_load_roundtrip_starts_a_clean_idle_runtime() {
    let mut app = EquilibriumApp::new();
    app.document.config.problem = EquilibriumProblemDraft::FixedPh {
        pressure_pa: "90000".into(),
        reference_pressure_pa: "101325".into(),
        target_enthalpy_j: "1234".into(),
        temperature_bounds: super::equilibrium_gui_model::PhTemperatureBoundsDraft::default(),
    };
    app.set_plot_series_visible("gas::old", false);
    let json = app.save_document_json().expect("document serializes");
    app.prepare_request()
        .expect("P,H request prepares without lookup");

    let mut restored = EquilibriumApp::new();
    restored
        .prepare_request()
        .expect("default request prepares before loading");
    restored
        .load_document_json(&json)
        .expect("versioned document loads");

    assert_eq!(restored.document, app.document);
    assert_eq!(
        restored.run_state(),
        super::equilibrium_gui_execution::EquilibriumGuiRunState::Idle
    );
    assert!(restored.prepared_request().is_none());
    assert!(restored.result_snapshot().is_none());
    assert!(restored.plot_series_visible("gas::old"));
}

#[test]
fn closing_and_reopening_the_window_preserves_only_the_editable_document() {
    let mut app = EquilibriumApp::new();
    app.document.config.problem = EquilibriumProblemDraft::FixedPt {
        pressure_pa: "90000".into(),
        reference_pressure_pa: "101325".into(),
        temperature: super::equilibrium_gui_model::TemperatureDraft::Point {
            temperature_k: "777".into(),
        },
    };
    let expected_document = app.document.clone();
    let context = egui::Context::default();
    let mut open = false;

    // Closing only hides the egui window. The app object remains the owner of
    // the document, while no accepted result or worker is manufactured by a
    // reopen operation.
    let _ = context.run_ui(Default::default(), |ctx| app.show(ctx, &mut open));
    assert!(!open);
    assert_eq!(app.document, expected_document);
    assert!(app.result_snapshot().is_none());

    open = true;
    let _ = context.run_ui(Default::default(), |ctx| app.show(ctx, &mut open));
    assert_eq!(app.document, expected_document);
    assert!(app.result_snapshot().is_none());
}

#[test]
fn equilibrium_tabs_expose_setup_and_advanced_groups_without_mixing_sections() {
    let app = Rc::new(RefCell::new(EquilibriumApp::new()));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });

    harness.run();
    harness.get_by_role_and_label(Role::Button, "Setup");
    harness.get_by_role_and_label(Role::Button, "Phase control");
    harness.get_by_role_and_label(Role::Button, "Libraries");
    harness.get_by_role_and_label(Role::Button, "Numerics");
    harness.get_by_role_and_label(Role::Button, "Output");
    harness.get_by_role_and_label(Role::Button, "Results");
    harness.get_by_label("P,T = const");
    harness.get_by_label("Explicit species");
    assert!(harness.query_by_label("Concrete backend").is_none());

    harness
        .get_by_role_and_label(Role::Button, "Phase control")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Phase policy");
    assert!(harness.query_by_label("Components and phases").is_none());

    harness
        .get_by_role_and_label(Role::Button, "Libraries")
        .click_accesskit();
    harness.run();
    harness.get_by_label("Library lookup");
    assert!(harness.query_by_label("Phase policy").is_none());
}

#[test]
fn candidate_preview_is_derived_and_becomes_stale_after_document_edit() {
    let mut app = EquilibriumApp::new();
    app.document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["C".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft {
            element_mode: GuiElementSearchMode::Exact,
            physical_states: Vec::new(),
            temperature_lower_k: "300".into(),
            temperature_upper_k: "1000".into(),
            max_candidates: "1".into(),
        },
        assignments: vec![super::equilibrium_gui_model::PhaseDraft::default()],
    };
    let repository = Arc::new(ThermoRepository::from_parts(
        vec![
            ("NASA_gas".into(), "CO".into()),
            ("NASA_gas".into(), "CO2".into()),
        ],
        HashMap::from([(
            "NASA_gas".into(),
            HashMap::from([
                ("CO".into(), json!({"T": [[200.0, 6000.0]]})),
                ("CO2".into(), json!({"T": [[200.0, 6000.0]]})),
            ]),
        )]),
        HashMap::from([
            (
                "C".into(),
                vec![
                    vec!["CO".into(), "NASA_gas".into()],
                    vec!["CO2".into(), "NASA_gas".into()],
                ],
            ),
            (
                "O".into(),
                vec![
                    vec!["CO".into(), "NASA_gas".into()],
                    vec!["CO2".into(), "NASA_gas".into()],
                ],
            ),
        ]),
        vec!["NASA_gas".into()],
        HashMap::new(),
        HashMap::new(),
        vec!["NASA_gas".into()],
        Vec::new(),
    ));

    app.preview_candidates(repository)
        .expect("candidate preview succeeds");
    assert_eq!(app.candidate_preview().unwrap().selected_count(), 1);
    let selected_key = app
        .candidate_preview()
        .unwrap()
        .rows()
        .iter()
        .find(|row| row.included())
        .and_then(|row| row.record_key())
        .expect("selected row has an exact record key")
        .to_string();
    assert!(app.assign_candidate_to_target_phase(&selected_key));
    if let EquilibriumInventoryDraft::ElementCandidates { assignments, .. } =
        &app.document.config.inventory
    {
        let assigned = assignments[0]
            .components
            .iter()
            .find(|component| component.substance == selected_key)
            .expect("selected candidate is present in the target phase");
        assert_eq!(assigned.source_library.as_deref(), Some("NASA_gas"));
    } else {
        panic!("candidate mode must retain phase assignments");
    }
    assert_eq!(app.candidate_initial_moles(&selected_key), Some("1.0"));
    assert!(app.set_candidate_initial_moles(&selected_key, "2.5"));
    assert_eq!(app.candidate_initial_moles(&selected_key), Some("2.5"));
    assert!(
        app.candidate_preview().is_some(),
        "dependent assignments must not invalidate the candidate query preview"
    );
    assert!(app.unassign_candidate(&selected_key));
    assert!(app.candidate_preview().is_some());
    if let EquilibriumInventoryDraft::ElementCandidates { elements, .. } =
        &mut app.document.config.inventory
    {
        elements[0] = "H".into();
    }
    assert!(app.candidate_preview().is_none());
}

#[test]
#[ignore = "requires the bundled local thermochemical catalog and element index"]
fn offline_local_element_candidate_story_keeps_catalogs_unchanged_and_renders_audit() {
    let before = local_gui_library_snapshot();
    let mut app = EquilibriumApp::new();
    app.document.config.inventory = EquilibriumInventoryDraft::ElementCandidates {
        elements: vec!["H".into(), "O".into()],
        candidate_policy: CandidatePolicyDraft {
            element_mode: GuiElementSearchMode::Exact,
            physical_states: vec![super::equilibrium_gui_model::GuiPhysicalState::Gas],
            temperature_lower_k: "300".into(),
            temperature_upper_k: "1000".into(),
            max_candidates: "20".into(),
        },
        assignments: vec![super::equilibrium_gui_model::PhaseDraft::default()],
    };
    app.document.config.lookup = super::equilibrium_gui_model::EquilibriumLookupDraft::Explicit {
        priority_libraries: vec!["NASA_gas".into()],
        permitted_libraries: vec!["NASA_gas".into()],
        explicit_search_instructions: Default::default(),
        search_in_nist: false,
    };

    app.preview_candidates(
        ThermoData::try_default_repository().expect("local thermochemical repository loads"),
    )
    .expect("real element candidate preview succeeds");

    let preview = app
        .candidate_preview()
        .expect("candidate preview remains current after discovery");
    assert_eq!(preview.requested_elements(), &["H", "O"]);
    assert!(preview.selected_count() > 0);
    assert!(preview.rows().iter().any(|row| row.included()));
    assert!(preview.rows().iter().all(|row| {
        !row.included() || (row.record_key().is_some() && !row.library().is_empty())
    }));
    let selected_keys = preview
        .rows()
        .iter()
        .filter(|row| row.included())
        .take(3)
        .filter_map(|row| row.record_key())
        .map(str::to_string)
        .collect::<Vec<_>>();
    assert!(
        selected_keys.len() >= 2,
        "the real catalog must expose at least two assignable H/O candidates"
    );
    // The preview is an input-discovery step, not a second solver model. Make
    // sure its selected record crosses the same assignment and request-build
    // boundary used by the eventual Run action.
    for selected_key in &selected_keys {
        if app.candidate_initial_moles(selected_key).is_none() {
            assert!(app.assign_candidate_to_target_phase(selected_key));
        }
    }
    app.prepare_request()
        .expect("assigned live candidate builds a canonical request");
    assert!(app.start_prepared_run());
    wait_for_gui_worker(&mut app);
    assert!(
        matches!(
            app.run_state(),
            super::equilibrium_gui_execution::EquilibriumGuiRunState::Completed { .. }
        ),
        "assigned element-mode request failed: state={:?}, error={:?}",
        app.run_state(),
        app.last_error()
    );
    assert!(app.result_snapshot().is_some());
    assert_eq!(before, local_gui_library_snapshot());

    let app = Rc::new(RefCell::new(app));
    let app_for_ui = Rc::clone(&app);
    let mut open = true;
    let mut harness = Harness::new_ui(move |ui| {
        app_for_ui.borrow_mut().show(ui.ctx(), &mut open);
    });
    harness.run();
    harness.get_by_label("Target phase");
    // Every included candidate has its own Assign button, so this label is
    // intentionally non-unique. The table header and target-phase control
    // are the stable accessibility anchors for this rendered audit.
    harness.get_by_label("Decision");
}
