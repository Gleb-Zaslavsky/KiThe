//! Helpers for rendering read-only GUI snapshots.
//!
//! These widgets present calculated reports and raw catalog payloads as
//! selectable/copyable text while keeping the underlying application state
//! immutable from direct user editing.

use eframe::egui;

/// Renders a multiline read-only snapshot that can still be selected and copied.
pub fn multiline(ui: &mut egui::Ui, text: &mut String, rows: usize) -> egui::Response {
    ui.add(
        egui::TextEdit::multiline(text)
            .desired_rows(rows)
            .desired_width(f32::INFINITY)
            .interactive(false)
            .code_editor(),
    )
}
