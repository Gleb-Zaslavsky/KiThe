use super::test_options::TestOptions;
use eframe::egui;

use std::ops::Range;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
pub struct HelpSectionNode {
    section_index: usize,
    children: Vec<usize>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct HelpDoc {
    pub(crate) raw_markdown: String,
    pub(crate) blocks: Vec<HelpBlock>,
    pub(crate) sections: Vec<HelpSection>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct HelpSection {
    pub(crate) title: String,
    pub(crate) level: usize,
    pub(crate) block_range: Range<usize>,
    pub(crate) search_blob: String,
}

#[derive(Clone, Debug, Default)]
pub(crate) enum HelpBlock {
    #[default]
    Empty,
    Heading {
        level: usize,
        text: String,
    },
    Paragraph(String),
    BulletList(Vec<String>),
    CodeBlock(String),
    Image {
        alt: String,
        path: String,
    },
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Default)]
pub(crate) enum HelpLanguage {
    #[default]
    Eng,
    Rus,
}

impl HelpDoc {
    pub fn hierarchical_sections(&self) -> Vec<HelpSectionNode> {
        let mut result = Vec::new();
        let mut current_parent_idx: Option<usize> = None;

        for (idx, section) in self.sections.iter().enumerate() {
            if section.level == 2 {
                result.push(HelpSectionNode {
                    section_index: idx,
                    children: Vec::new(),
                });
                current_parent_idx = Some(result.len() - 1);
            } else if section.level == 1 {
                if let Some(parent_idx) = current_parent_idx {
                    result[parent_idx].children.push(idx);
                } else {
                    result.push(HelpSectionNode {
                        section_index: idx,
                        children: Vec::new(),
                    });
                    current_parent_idx = Some(result.len() - 1);
                }
            } else {
                // Levels >2 are not considered for top-level hierarchy in current behavior.
                // They can be handled similarly if needed in future.
            }
        }

        result
    }
}

impl HelpLanguage {
    fn file_name(self) -> &'static str {
        match self {
            Self::Eng => "help_eng.md",
            Self::Rus => "help_rus.md",
        }
    }
}

impl TestOptions {
    pub fn show_help_window_ui(&mut self, ctx: &egui::Context) {
        if !self.show_help_window {
            return;
        }

        egui_extras::install_image_loaders(ctx);

        let mut open = self.show_help_window;
        egui::Window::new("Help")
            .open(&mut open)
            .resizable(true)
            .default_width(980.0)
            .default_height(700.0)
            .show(ctx, |ui| {
                self.render_help_toolbar(ui);

                if let Some(err) = &self.help_error {
                    ui.colored_label(egui::Color32::RED, err);
                }
                if let Some(status) = &self.help_status {
                    ui.colored_label(egui::Color32::LIGHT_GREEN, status);
                }

                ui.separator();

                if self.help_edit_mode {
                    self.render_help_editor(ui);
                } else {
                    self.render_help_browser(ui);
                }
            });
        self.show_help_window = open;
    }

    fn render_help_toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.label("Language:");

            if ui
                .selectable_label(self.help_lang == HelpLanguage::Eng, "ENG")
                .clicked()
            {
                self.help_lang = HelpLanguage::Eng;
                let _ = self.reload_help();
            }
            if ui
                .selectable_label(self.help_lang == HelpLanguage::Rus, "RUS")
                .clicked()
            {
                self.help_lang = HelpLanguage::Rus;
                let _ = self.reload_help();
            }

            ui.separator();
            ui.label("Search:");
            let search = ui.add_sized(
                [240.0, 24.0],
                egui::TextEdit::singleline(&mut self.help_search_query)
                    .hint_text("section, command, term..."),
            );
            if search.changed() && self.help_selected_section.is_none() {
                self.help_selected_section = self.filtered_section_indices().first().copied();
            }

            ui.separator();
            if ui.button("Reload").clicked() {
                let _ = self.reload_help();
            }

            if !self.help_edit_mode {
                if ui.button("Edit").clicked() {
                    self.help_edit_mode = true;
                    if self.help_editor_buffer.is_empty() {
                        self.help_editor_buffer = self.parsed_help.raw_markdown.clone();
                    }
                    self.help_status = Some("Edit mode enabled".to_string());
                }
            } else {
                if ui.button("Save").clicked() {
                    if let Err(err) = self.save_help() {
                        self.help_error = Some(err);
                    }
                }
                if ui.button("Cancel").clicked() {
                    self.help_edit_mode = false;
                    self.help_editor_buffer = self.parsed_help.raw_markdown.clone();
                    self.help_status = Some("Changes discarded".to_string());
                }
            }

            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.button("Close").clicked() {
                    self.show_help_window = false;
                }
            });
        });
    }

    fn render_help_editor(&mut self, ui: &mut egui::Ui) {
        ui.label(format!("Editing: {}", self.current_help_path().display()));
        ui.add_space(6.0);
        egui::ScrollArea::vertical()
            .id_salt("help_editor_scroll")
            .show(ui, |ui| {
                ui.add(
                    egui::TextEdit::multiline(&mut self.help_editor_buffer)
                        .desired_rows(32)
                        .desired_width(f32::INFINITY)
                        .code_editor(),
                );
            });
    }

    fn render_help_browser(&mut self, ui: &mut egui::Ui) {
        let filtered_sections = self.filtered_section_indices();
        if filtered_sections.is_empty() {
            ui.label("No sections match the current search query.");
            return;
        }

        let selected = self
            .help_selected_section
            .filter(|idx| filtered_sections.contains(idx))
            .unwrap_or(filtered_sections[0]);
        self.help_selected_section = Some(selected);

        ui.columns(2, |columns| {
            columns[0].set_min_width(240.0);
            columns[0].heading("Sections");
            columns[0].separator();
            egui::ScrollArea::vertical()
                .id_salt("help_sections_list")
                .show(&mut columns[0], |ui| {
                    if self.help_search_query.trim().is_empty() {
                        for node in self.parsed_help.hierarchical_sections().iter() {
                            let section = &self.parsed_help.sections[node.section_index];
                            let is_open = self.help_expanded_sections.contains(&node.section_index);

                            ui.horizontal(|ui| {
                                if !node.children.is_empty() {
                                    if ui.small_button(if is_open { "▾" } else { "▸" }).clicked()
                                    {
                                        if is_open {
                                            self.help_expanded_sections.remove(&node.section_index);
                                        } else {
                                            self.help_expanded_sections.insert(node.section_index);
                                        }
                                    }
                                } else {
                                    ui.add_space(14.0);
                                }

                                if ui
                                    .selectable_label(
                                        self.help_selected_section == Some(node.section_index),
                                        &section.title,
                                    )
                                    .clicked()
                                {
                                    self.help_selected_section = Some(node.section_index);
                                }
                            });

                            if is_open {
                                for &child_idx in node.children.iter() {
                                    let child = &self.parsed_help.sections[child_idx];
                                    ui.horizontal(|ui| {
                                        ui.add_space(18.0);
                                        if ui
                                            .selectable_label(
                                                self.help_selected_section == Some(child_idx),
                                                &child.title,
                                            )
                                            .clicked()
                                        {
                                            self.help_selected_section = Some(child_idx);
                                        }
                                    });
                                }
                            }
                        }
                    } else {
                        for idx in filtered_sections.iter().copied() {
                            let section = &self.parsed_help.sections[idx];
                            let indent = if section.level > 1 { "  " } else { "" };
                            let label = format!("{indent}{}", section.title);
                            if ui
                                .selectable_label(self.help_selected_section == Some(idx), label)
                                .clicked()
                            {
                                self.help_selected_section = Some(idx);
                            }
                        }
                    }
                });

            columns[1].heading(if self.help_search_query.trim().is_empty() {
                self.parsed_help.sections[selected].title.clone()
            } else {
                format!("Search results: {}", filtered_sections.len())
            });
            columns[1].separator();
            egui::ScrollArea::vertical()
                .id_salt("help_section_content")
                .show(&mut columns[1], |ui| {
                    if self.help_search_query.trim().is_empty() {
                        self.render_section(ui, selected, false);
                    } else {
                        for idx in filtered_sections.iter().copied() {
                            self.render_section(ui, idx, true);
                            ui.add_space(12.0);
                            ui.separator();
                        }
                    }
                });
        });
    }

    fn render_section(&self, ui: &mut egui::Ui, section_index: usize, search_mode: bool) {
        let section = &self.parsed_help.sections[section_index];
        if search_mode {
            ui.label(egui::RichText::new(&section.title).strong().size(24.0));
            ui.add_space(4.0);
        }

        for block in &self.parsed_help.blocks[section.block_range.clone()] {
            if search_mode && !block_matches_query(block, &self.help_search_query) {
                continue;
            }
            render_help_block(ui, block, self.current_help_path().parent());
        }
    }

    pub(crate) fn current_help_path(&self) -> PathBuf {
        PathBuf::from("src/assets").join(self.help_lang.file_name())
    }

    pub(crate) fn reload_help(&mut self) -> Result<(), String> {
        let path = self.current_help_path();
        match std::fs::read_to_string(&path) {
            Ok(contents) => {
                self.parsed_help = parse_help_markdown(&contents);
                self.help_editor_buffer = contents;
                self.help_error = None;
                self.help_status = Some(format!("Loaded {}", path.display()));
                self.help_selected_section = self.filtered_section_indices().first().copied();
                Ok(())
            }
            Err(err) => {
                self.parsed_help = HelpDoc::default();
                self.help_error = Some(format!("Failed to read {}: {}", path.display(), err));
                Err(err.to_string())
            }
        }
    }

    fn save_help(&mut self) -> Result<(), String> {
        let path = self.current_help_path();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|err| err.to_string())?;
        }
        std::fs::write(&path, &self.help_editor_buffer).map_err(|err| err.to_string())?;
        self.help_edit_mode = false;
        self.help_status = Some(format!("Saved {}", path.display()));
        self.reload_help()
    }

    fn filtered_section_indices(&self) -> Vec<usize> {
        let query = self.help_search_query.trim();
        if query.is_empty() {
            return (0..self.parsed_help.sections.len()).collect();
        }

        self.parsed_help
            .sections
            .iter()
            .enumerate()
            .filter_map(|(idx, section)| {
                if contains_case_insensitive(&section.search_blob, query) {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect()
    }
}

fn parse_help_markdown(src: &str) -> HelpDoc {
    let mut blocks = Vec::new();
    let mut paragraph_lines: Vec<String> = Vec::new();
    let mut list_items: Vec<String> = Vec::new();
    let mut code_lines: Vec<String> = Vec::new();
    let mut in_code_block = false;

    let flush_paragraph = |blocks: &mut Vec<HelpBlock>, paragraph_lines: &mut Vec<String>| {
        if !paragraph_lines.is_empty() {
            blocks.push(HelpBlock::Paragraph(paragraph_lines.join("\n")));
            paragraph_lines.clear();
        }
    };
    let flush_list = |blocks: &mut Vec<HelpBlock>, list_items: &mut Vec<String>| {
        if !list_items.is_empty() {
            blocks.push(HelpBlock::BulletList(std::mem::take(list_items)));
        }
    };
    let flush_code = |blocks: &mut Vec<HelpBlock>, code_lines: &mut Vec<String>| {
        if !code_lines.is_empty() {
            blocks.push(HelpBlock::CodeBlock(code_lines.join("\n")));
            code_lines.clear();
        }
    };

    for raw_line in src.lines() {
        let line = raw_line.trim_end_matches('\r');
        let trimmed = line.trim();

        if in_code_block {
            if trimmed.starts_with("```") {
                flush_code(&mut blocks, &mut code_lines);
                in_code_block = false;
            } else {
                code_lines.push(line.to_string());
            }
            continue;
        }

        if trimmed.starts_with("```") {
            flush_paragraph(&mut blocks, &mut paragraph_lines);
            flush_list(&mut blocks, &mut list_items);
            in_code_block = true;
            continue;
        }

        if trimmed.is_empty() {
            flush_paragraph(&mut blocks, &mut paragraph_lines);
            flush_list(&mut blocks, &mut list_items);
            continue;
        }

        if let Some((level, text)) = parse_heading(trimmed) {
            flush_paragraph(&mut blocks, &mut paragraph_lines);
            flush_list(&mut blocks, &mut list_items);
            blocks.push(HelpBlock::Heading { level, text });
            continue;
        }

        if let Some((alt, path)) = parse_image(trimmed) {
            flush_paragraph(&mut blocks, &mut paragraph_lines);
            flush_list(&mut blocks, &mut list_items);
            blocks.push(HelpBlock::Image { alt, path });
            continue;
        }

        if let Some(item) = parse_list_item(trimmed) {
            flush_paragraph(&mut blocks, &mut paragraph_lines);
            list_items.push(item);
            continue;
        }

        flush_list(&mut blocks, &mut list_items);
        paragraph_lines.push(line.to_string());
    }

    flush_paragraph(&mut blocks, &mut paragraph_lines);
    flush_list(&mut blocks, &mut list_items);
    flush_code(&mut blocks, &mut code_lines);

    let mut sections = Vec::new();
    let headings: Vec<(usize, usize, String)> = blocks
        .iter()
        .enumerate()
        .filter_map(|(idx, block)| match block {
            HelpBlock::Heading { level, text } if *level <= 2 => Some((idx, *level, text.clone())),
            _ => None,
        })
        .collect();

    if headings.is_empty() {
        sections.push(HelpSection {
            title: "Help".to_string(),
            level: 1,
            block_range: 0..blocks.len(),
            search_blob: blocks
                .iter()
                .map(block_search_text)
                .collect::<Vec<_>>()
                .join("\n"),
        });
    } else {
        for (position, (start_idx, level, title)) in headings.iter().enumerate() {
            let end_idx = headings
                .get(position + 1)
                .map(|(next_idx, _, _)| *next_idx)
                .unwrap_or(blocks.len());
            let blob = blocks[*start_idx..end_idx]
                .iter()
                .map(block_search_text)
                .collect::<Vec<_>>()
                .join("\n");
            sections.push(HelpSection {
                title: title.clone(),
                level: *level,
                block_range: *start_idx..end_idx,
                search_blob: blob,
            });
        }
    }

    HelpDoc {
        raw_markdown: src.to_string(),
        blocks,
        sections,
    }
}

fn parse_heading(line: &str) -> Option<(usize, String)> {
    let level = line.chars().take_while(|ch| *ch == '#').count();
    if !(1..=6).contains(&level) {
        return None;
    }
    let text = line[level..].trim();
    if text.is_empty() {
        None
    } else {
        Some((level, text.to_string()))
    }
}

fn parse_image(line: &str) -> Option<(String, String)> {
    if !line.starts_with("![") {
        return None;
    }
    let alt_end = line.find("](")?;
    let path_end = line.rfind(')')?;
    if path_end <= alt_end + 2 {
        return None;
    }
    Some((
        line[2..alt_end].trim().to_string(),
        line[alt_end + 2..path_end].trim().to_string(),
    ))
}

fn parse_list_item(line: &str) -> Option<String> {
    if let Some(item) = line.strip_prefix("- ") {
        return Some(item.trim().to_string());
    }
    if let Some(item) = line.strip_prefix("* ") {
        return Some(item.trim().to_string());
    }
    let mut chars = line.chars();
    let mut saw_digit = false;
    while let Some(ch) = chars.next() {
        if ch.is_ascii_digit() {
            saw_digit = true;
            continue;
        }
        if ch == '.' && saw_digit {
            let rest = chars.as_str().trim();
            if !rest.is_empty() {
                return Some(rest.to_string());
            }
        }
        break;
    }
    None
}

fn render_help_block(ui: &mut egui::Ui, block: &HelpBlock, base_dir: Option<&Path>) {
    match block {
        HelpBlock::Empty => {}
        HelpBlock::Heading { level, text } => {
            let size = match level {
                1 => 28.0,
                2 => 22.0,
                _ => 18.0,
            };
            ui.add_space(4.0);
            ui.label(egui::RichText::new(text).strong().size(size));
            ui.add_space(2.0);
        }
        HelpBlock::Paragraph(text) => {
            ui.label(text);
            ui.add_space(6.0);
        }
        HelpBlock::BulletList(items) => {
            for item in items {
                ui.horizontal_wrapped(|ui| {
                    ui.label("?");
                    ui.label(item);
                });
            }
            ui.add_space(6.0);
        }
        HelpBlock::CodeBlock(code) => {
            egui::Frame::group(ui.style()).show(ui, |ui| {
                ui.label(egui::RichText::new(code).monospace());
            });
            ui.add_space(6.0);
        }
        HelpBlock::Image { alt, path } => {
            let resolved = resolve_image_path(base_dir, path);
            if resolved.exists() {
                let uri = format!("file://{}", resolved.to_string_lossy().replace('\\', "/"));
                ui.add(
                    egui::Image::from_uri(uri)
                        .alt_text(alt)
                        .shrink_to_fit()
                        .max_height(420.0),
                );
                if !alt.is_empty() {
                    ui.small(alt);
                }
            } else {
                ui.colored_label(
                    egui::Color32::YELLOW,
                    format!("Image not found: {}", resolved.display()),
                );
            }
            ui.add_space(8.0);
        }
    }
}

fn resolve_image_path(base_dir: Option<&Path>, path: &str) -> PathBuf {
    let candidate = PathBuf::from(path);
    if candidate.is_absolute() {
        candidate
    } else if let Some(base_dir) = base_dir {
        base_dir.join(candidate)
    } else {
        candidate
    }
}

fn block_search_text(block: &HelpBlock) -> String {
    match block {
        HelpBlock::Empty => String::new(),
        HelpBlock::Heading { text, .. } => text.clone(),
        HelpBlock::Paragraph(text) => text.clone(),
        HelpBlock::BulletList(items) => items.join("\n"),
        HelpBlock::CodeBlock(code) => code.clone(),
        HelpBlock::Image { alt, path } => format!("{alt}\n{path}"),
    }
}

fn block_matches_query(block: &HelpBlock, query: &str) -> bool {
    let query = query.trim();
    query.is_empty() || contains_case_insensitive(&block_search_text(block), query)
}

fn contains_case_insensitive(haystack: &str, needle: &str) -> bool {
    haystack.to_lowercase().contains(&needle.to_lowercase())
}
