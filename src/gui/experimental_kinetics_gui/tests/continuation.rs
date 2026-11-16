use super::support::{empty_app, import_csv, with_app};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
use egui_kittest::kittest::Queryable;
use tempfile::tempdir;

#[test]
fn story_continue_experiment_from_second_file_and_reject_bad_time_order_atomically() {
    let dir = tempdir().unwrap();
    let first = dir.path().join("segment_a.csv");
    let next = dir.path().join("segment_b.csv");
    let invalid = dir.path().join("segment_invalid.csv");
    std::fs::write(
        &first,
        "t,T,m\n0.0,300.0,10.0\n1.0,301.0,9.0\n2.0,302.0,8.0\n",
    )
    .unwrap();
    // A reset clock is shifted to the next target sample time by the append workflow.
    std::fs::write(
        &next,
        "t,T,m\n0.0,303.0,7.0\n1.0,304.0,6.0\n2.0,305.0,5.0\n",
    )
    .unwrap();
    std::fs::write(
        &invalid,
        "t,T,m\n0.0,303.0,7.0\n2.0,304.0,6.0\n1.0,305.0,5.0\n",
    )
    .unwrap();

    with_app(empty_app(), |harness| {
        import_csv(harness, &first, "segment_a", None, None);
        let before = super::support::read_column(&harness.state().app.model, "segment_a", "t");
        assert_eq!(before, vec![0.0, 1.0, 2.0]);

        super::support::click_menu_entry(harness, "File Manager", "Append to experiment");
        harness.get_by_label("[F] segment_invalid.csv").click();
        harness.run();
        harness.get_by_label("Preview continuation").click();
        harness.run();
        assert!(
            harness
                .state()
                .app
                .new_experiment_dialog
                .last_error
                .is_some()
        );
        assert_eq!(harness.state().app.model.series.len(), 1);
        assert_eq!(
            super::support::read_column(&harness.state().app.model, "segment_a", "t"),
            before,
            "failed preview must not modify its source experiment"
        );

        harness.get_by_label("[F] segment_b.csv").click();
        harness.run();
        harness.get_by_label("Preview continuation").click();
        harness.run();
        assert!(
            harness
                .state()
                .app
                .new_experiment_dialog
                .last_error
                .is_none()
        );
        harness.get_by_label("Add combined experiment").click();
        harness.run();
        harness.run();

        let state = harness.state();
        let model = &state.app.model;
        assert_eq!(model.series.len(), 2, "append keeps the original segment");
        let combined_id = "segment_a_continued";
        let times = super::support::read_column(model, combined_id, "t");
        assert_eq!(times, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        assert_eq!(
            super::support::read_column(model, combined_id, "T"),
            vec![300.0, 301.0, 302.0, 303.0, 304.0, 305.0]
        );
        assert_eq!(
            super::support::read_column(model, combined_id, "m"),
            vec![10.0, 9.0, 8.0, 7.0, 6.0, 5.0]
        );
        let dataset = &model
            .series
            .get_experiment_by_id(combined_id)
            .unwrap()
            .dataset;
        assert_eq!(dataset.schema.columns["t"].unit, Unit::Second);
        assert_eq!(dataset.schema.time.as_deref(), Some("t"));
        assert_eq!(dataset.schema.temperature.as_deref(), Some("T"));
        assert_eq!(dataset.schema.mass.as_deref(), Some("m"));
        assert!(
            model
                .plots
                .iter()
                .any(|curve| curve.experiment_id == combined_id)
        );
        println!(
            "[STORY][continuation] source_rows=3 appended_rows=3 inferred_dt=1s shifted_boundary=3s combined_rows={} columns=t,T,m invalid_out_of_order_preview=rejected source_preserved=true plots_refreshed=true",
            times.len()
        );
    });
}
