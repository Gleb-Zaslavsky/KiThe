//! Static contextual help for the equilibrium calculator.
//!
//! The lookup is deliberately pure and keyed by stable identifiers. It does
//! not inspect the editable document, repository, or solver state.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumHelpLanguage {
    English,
    Russian,
}

/// A stable tooltip identifier registered by the equilibrium help catalog.
///
/// The ID is intentionally still a string at the resource/UI boundary: egui
/// labels, Markdown help, and future translated resources need a durable
/// external name. The descriptor catalog is the single source for registered
/// control IDs. Actual rendered-widget coverage belongs to story tests, not
/// to a manually maintained second inventory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EquilibriumTooltipDescriptor {
    pub key: EquilibriumHelpKey,
}

/// Validated tooltip key used by internal help consumers.
///
/// The constructor checks the canonical catalog, so an arbitrary string can
/// enter only at the UI/resource boundary. Keeping the ID inside this type
/// makes it possible to migrate call sites incrementally without duplicating
/// the catalog as a large enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EquilibriumHelpKey(&'static str);

impl EquilibriumHelpKey {
    /// Creates a key only when it is registered in the canonical catalog.
    pub fn try_from(id: &str) -> Option<Self> {
        tooltip_descriptor(id).map(|descriptor| descriptor.key)
    }

    pub const fn as_str(self) -> &'static str {
        self.0
    }
}

macro_rules! tooltip_catalog {
    ($($id:literal),+ $(,)?) => {
        &[
            $(EquilibriumTooltipDescriptor { key: EquilibriumHelpKey($id) }),+
        ]
    };
}

/// Canonical tooltip catalog. It is the only declaration of tooltip identity;
/// tests derive catalog completeness and duplicate checks from this table.
pub const EQUILIBRIUM_TOOLTIP_CATALOG: &[EquilibriumTooltipDescriptor] = tooltip_catalog![
    "tab.setup",
    "tab.phase_control",
    "tab.libraries",
    "tab.numerics",
    "tab.output",
    "tab.results",
    "help.language",
    "action.validate",
    "action.prepare",
    "action.run",
    "action.cancel",
    "action.embedded_plot",
    "action.kithe_plot",
    "mode.pt",
    "mode.ph",
    "inventory.explicit",
    "inventory.elements",
    "inventory.state",
    "inventory.model",
    "candidate.matching",
    "candidate.policy",
    "candidate.states",
    "candidate.bounds",
    "candidate.max",
    "candidate.preview",
    "candidate.cancel",
    "candidate.target",
    "lookup.priority",
    "lookup.permitted",
    "lookup.nist",
    "solver.budget",
    "solver.scaling",
    "solver.trace_override",
    "solver.trace_seed",
    "output.basis",
    "output.plot_target",
    "output.pchip",
    "output.clamp",
    "output.interpolation",
    "output.scale",
    "temperature.point",
    "temperature.range",
    "phase.initial_set",
    "phase.fixed",
    "phase.bounded",
    "solver.ph_route",
    "diagnostics.timing",
    "diagnostics.phase_trace",
    "diagnostics.range_trace",
    "diagnostics.keq",
    "list.up",
    "list.down",
    "list.remove",
    "list.add_phase",
    "list.add_component",
    "list.assign",
    "list.unassign",
    "list.add_backend",
    "list.add_library",
    "preset.ideal_gas",
    "output.visible_series",
    "field.pressure",
    "field.reference_pressure",
    "field.target_enthalpy",
    "field.phase_id",
    "field.substance",
    "field.moles",
    "field.element",
    "field.phase_epsilon",
    "field.dg_create",
    "field.dg_keep",
    "field.phase_iterations",
    "field.tolerance",
    "field.max_iterations",
    "field.cascade_attempts",
    "field.cascade_iterations",
    "field.cascade_total",
    "field.trace_floor",
    "field.trace_fraction",
    "field.trace_minimum_floor",
    "field.display_points",
    "field.display_cutoff",
    "lookup.priority_entry",
    "lookup.permitted_entry",
    "state.gas",
    "state.liquid",
    "state.solid",
    "model.ideal_gas",
    "model.ideal_solution",
    "model.pure_condensed",
    "matching.subset",
    "matching.exact",
    "basis.moles",
    "basis.fractions",
    "basis.phase_totals",
    "scale.linear",
    "scale.log",
    "interpolation.linear",
    "interpolation.log",
    "language.english",
    "language.russian",
    "list.add_element",
    "candidate.state.gas",
    "candidate.state.liquid",
    "candidate.state.solid",
    "candidate.state.condensed",
    "lookup.catalog",
    "lookup.load_catalog",
    "lookup.default",
    "lookup.explicit",
    "phase.initial_positive",
    "phase.initial_all",
    "solver.ph_auto",
    "solver.ph_monolithic",
    "solver.ph_nested",
    "solver.production",
    "solver.backend",
    "solver.custom",
    "solver.trace_absolute",
    "solver.trace_relative",
    "diagnostics.phase_off",
    "diagnostics.phase_summary",
    "diagnostics.phase_lifecycle",
    "diagnostics.phase_detailed",
    "diagnostics.range_endpoints",
    "diagnostics.range_transitions",
    "diagnostics.range_every_point",
    "diagnostics.keq_off",
    "diagnostics.keq_when_applicable",
    "diagnostics.keq_required",
    "output.plot_none",
    "output.plot_embedded",
    "output.plot_kithe_plot",
    "output.plot_both",
    "output.table_density",
    "output.table_detailed",
    "output.table_compact",
];

/// Registered keys that intentionally retain English short text in the
/// Russian interface while their resource translation is still pending.
///
/// Keep this list small and explicit. The catalog test rejects every other
/// English/Russian text collision, so a newly added control cannot silently
/// lose localization.
pub const RUSSIAN_TOOLTIP_FALLBACK_KEYS: &[&str] = &[];

/// Finds a descriptor by its stable resource/UI identifier.
pub fn tooltip_descriptor(control_key: &str) -> Option<&'static EquilibriumTooltipDescriptor> {
    EQUILIBRIUM_TOOLTIP_CATALOG
        .iter()
        .find(|descriptor| descriptor.key.as_str() == control_key)
}

/// Typed counterpart of [`tooltip`].
pub fn tooltip_typed(language: EquilibriumHelpLanguage, key: EquilibriumHelpKey) -> &'static str {
    tooltip(language, key.as_str())
}

/// Iterates stable IDs without exposing a second inventory that can drift from
/// the canonical descriptor catalog.
pub fn registered_tooltip_ids() -> impl Iterator<Item = &'static str> {
    EQUILIBRIUM_TOOLTIP_CATALOG
        .iter()
        .map(|descriptor| descriptor.key.as_str())
}

const HELP_UNAVAILABLE: &str = "Help text is unavailable for this tab.";

/// Extracts one top-level Help tab from a Markdown resource.
///
/// Keeping parsing separate lets tests exercise a missing localized section
/// without mutating compile-time resources or constructing a GUI worker.
fn section_from_document<'a>(document: &'a str, tab_key: &str) -> Option<&'a str> {
    let marker = format!("## {tab_key}\n");
    document
        .split_once(&marker)
        .and_then(|(_, rest)| rest.split_once("\n## ").map(|(section, _)| section))
        .or_else(|| document.split_once(&marker).map(|(_, rest)| rest))
}

fn document_for(language: EquilibriumHelpLanguage) -> &'static str {
    match language {
        EquilibriumHelpLanguage::English => include_str!("../assets/equilibrium_help_eng.md"),
        EquilibriumHelpLanguage::Russian => include_str!("../assets/equilibrium_help_rus.md"),
    }
}

pub fn text(language: EquilibriumHelpLanguage, tab_key: &str) -> &'static str {
    section_from_document(document_for(language), tab_key).unwrap_or(HELP_UNAVAILABLE)
}

/// Returns localized help and falls back to English when a section is absent.
pub fn text_with_fallback(language: EquilibriumHelpLanguage, tab_key: &str) -> &'static str {
    section_from_document(document_for(language), tab_key)
        .or_else(|| section_from_document(document_for(EquilibriumHelpLanguage::English), tab_key))
        .unwrap_or(HELP_UNAVAILABLE)
}

/// Reads the structured subset of short tooltips from the paired resources.
/// The parser deliberately mirrors the full-help format and returns borrowed
/// `include_str!` slices, so lookup stays offline, allocation-free, and static.
fn resource_tooltip(language: EquilibriumHelpLanguage, control_key: &str) -> Option<&'static str> {
    let document = match language {
        EquilibriumHelpLanguage::English => {
            include_str!("../assets/equilibrium_tooltips_eng.md")
        }
        EquilibriumHelpLanguage::Russian => {
            include_str!("../assets/equilibrium_tooltips_rus.md")
        }
    };
    let marker = format!("## {control_key}\n");
    document
        .split_once(&marker)
        .map(|(_, rest)| rest.split_once("\n## ").map_or(rest, |(text, _)| text))
        .map(str::trim)
        .filter(|text| !text.is_empty())
}

/// Returns the short hover explanation for a stable control identifier.
pub fn tooltip(language: EquilibriumHelpLanguage, control_key: &str) -> &'static str {
    resource_tooltip(language, control_key).unwrap_or("Equilibrium calculator control.")
}

/// Compatibility-only lookup retained while external callers migrate to the
/// structured tooltip resources. Internal GUI rendering never uses this path.
#[deprecated(note = "use tooltip; legacy Rust match will be removed")]
pub fn legacy_tooltip(language: EquilibriumHelpLanguage, control_key: &str) -> &'static str {
    if let Some(resource_text) = resource_tooltip(language, control_key) {
        return resource_text;
    }

    let english = match control_key {
        "field.pressure" => {
            "System pressure in pascals; changing it invalidates the prepared request."
        }
        "field.reference_pressure" => {
            "Gas activity reference pressure; it is distinct from system pressure."
        }
        "field.phase_id" => {
            "Stable identifier used to distinguish this phase in reports and plots."
        }
        "field.substance" => {
            "Thermochemical component identifier resolved through the library policy."
        }
        "field.display_cutoff" => {
            "Presentation cutoff for small mole fractions; it does not change solved data."
        }
        "field.generic" => "Equilibrium calculator control.",
        "mode.pt" => "Solve at fixed pressure and temperature, or over a temperature range.",
        "mode.ph" => {
            "Solve at fixed pressure and total enthalpy; temperature is determined by the solver."
        }
        "inventory.explicit" => "Enter the phase-qualified substances directly.",
        "inventory.elements" => {
            "Search selected catalogs for candidates matching the chosen elements."
        }
        "preset.ideal_gas" => "Replace the inventory with a small valid ideal-gas example.",
        "diagnostics.timing" => "Retain stage and per-point timing in the diagnostic report.",
        "diagnostics.phase_trace" => "Select how much phase lifecycle evidence is retained.",
        "diagnostics.range_trace" => {
            "Select which temperature-range lifecycle points are retained."
        }
        "diagnostics.keq" => "Use the independent equilibrium-constant validator when applicable.",
        "phase.initial_set" => {
            "Choose whether the initial active phase set follows inventory or all declared candidates."
        }
        "solver.ph_route" => "Choose Auto, monolithic, or nested temperature handling for P,H.",
        "output.resampling" => "Enable display-only PCHIP resampling for temperature ranges.",
        "candidate.matching" => {
            "Choose subset matching or require exactly the requested element set."
        }
        "candidate.bounds" => {
            "Limit candidate records to this thermochemistry temperature interval."
        }
        "candidate.preview" => {
            "Resolve and display candidates without starting equilibrium solving."
        }
        _ => "Equilibrium calculator control.",
    };
    if language == EquilibriumHelpLanguage::English {
        return english;
    }
    match control_key {
        "tab.setup" => "Задайте задачу, состав, фазы и начальные количества.",
        "tab.phase_control" => "Настройте активацию, удаление фаз, TPD и гистерезис.",
        "tab.libraries" => "Выберите источники данных и политику offline/fallback.",
        "tab.numerics" => "Настройте каскад нелинейных солверов и численную защиту.",
        "tab.output" => "Настройте диагностику, представление результатов и графики.",
        "tab.results" => "Просмотрите неизменяемый принятый snapshot и его доказательства.",
        "help.language" => "Выберите язык справки и всплывающих подсказок.",
        "action.validate" => "Проверьте документ без запуска расчета.",
        "action.prepare" => "Разрешите данные и соберите канонический запрос.",
        "action.run" => "Запустите подготовленный запрос в фоновом worker-е.",
        "action.cancel" => "Отмените расчет; частичный результат не публикуется.",
        "action.embedded_plot" => "Откройте график принятого snapshot без нового расчета.",
        "action.kithe_plot" => "Откройте KiThePlot для принятого snapshot.",
        "field.pressure" => {
            "Давление системы в паскалях; изменение отменяет подготовленный запрос."
        }
        "field.reference_pressure" => {
            "Опорное давление активности газа; оно отличается от давления системы."
        }
        "field.target_enthalpy" => "Целевая полная энтальпия в джоулях для ограничения P,H.",
        "field.phase_id" => "Стабильный идентификатор фазы в отчетах и графиках.",
        "field.substance" => "Идентификатор вещества, разрешаемый политикой библиотек.",
        "field.moles" => {
            "Начальное количество в молях; оно должно быть конечным и неотрицательным."
        }
        "field.phase_epsilon" => "Малый порог количества фазы для решений жизненного цикла.",
        "field.dg_create" => "Порог TPD, ниже которого отсутствующая фаза может появиться.",
        "field.dg_keep" => "Порог TPD, управляющий сохранением активной фазы.",
        "field.phase_iterations" => "Максимальное число итераций жизненного цикла фаз.",
        "field.tolerance" => "Необязательная замена допуска приемки нелинейного солвера.",
        "field.max_iterations" => "Необязательное ограничение итераций backend-а.",
        "field.display_points" => "Число точек только для отображения диапазона P,T.",
        "field.display_cutoff" => "Порог отображения малых долей; решенные данные не меняются.",
        "field.generic" => "Элемент управления калькулятора равновесия.",
        "output.pchip" => "Включите PCHIP только для отображения диапазона P,T.",
        "output.clamp" => "Ограничьте интерполированные значения решенным диапазоном.",
        "temperature.point" => "Решить одну точку при фиксированной температуре в кельвинах.",
        "temperature.range" => "Решить упорядоченную температурную сетку в кельвинах.",
        "field.element" => "Химический символ элемента для поиска кандидатов.",
        "field.trace_floor" => "Абсолютный положительный seed следового вещества в молях.",
        "field.trace_fraction" => "Следовой seed как доля наибольшего начального количества.",
        "field.trace_minimum_floor" => "Нижняя граница относительного следового seed.",
        "lookup.priority_entry" => "Имя библиотеки в упорядоченном списке приоритетов поиска.",
        "lookup.permitted_entry" => "Имя библиотеки из замкнутого набора разрешенных источников.",
        _ => english,
    }
}

/// Returns a tooltip only when the stable key is registered and resolves to
/// specific text. Production rendering uses [`tooltip`] so an incomplete
/// catalog never blocks the calculator; tests use this stricter contract to
/// prevent the generic fallback from hiding a missing entry.
pub fn tooltip_strict(
    language: EquilibriumHelpLanguage,
    control_key: &str,
) -> Option<&'static str> {
    let key = EquilibriumHelpKey::try_from(control_key)?;
    let resolved = tooltip_typed(language, key);
    if resolved != "Equilibrium calculator control." {
        return Some(resolved);
    }
    let option = option_tooltip(language, control_key);
    (option != "Equilibrium selection.").then_some(option)
}

/// Maps shared editor labels to the short explanations used by text fields.
pub fn field_tooltip(language: EquilibriumHelpLanguage, label: &str) -> &'static str {
    let key = match label {
        "Pressure [Pa]" => "field.pressure",
        "Reference pressure [Pa]" => "field.reference_pressure",
        "Target total enthalpy [J]" => "field.target_enthalpy",
        "Phase id" => "field.phase_id",
        "Substance" => "field.substance",
        "Moles" => "field.moles",
        "Element" => "field.element",
        "Phase epsilon" => "field.phase_epsilon",
        "Creation driving force" => "field.dg_create",
        "Keep driving force" => "field.dg_keep",
        "Maximum phase iterations" => "field.phase_iterations",
        "Tolerance override" => "field.tolerance",
        "Max iterations override" => "field.max_iterations",
        "Cascade max attempts" => "field.cascade_attempts",
        "Cascade iterations per attempt" => "field.cascade_iterations",
        "Cascade total iterations" => "field.cascade_total",
        "Trace floor" => "field.trace_floor",
        "Trace fraction" => "field.trace_fraction",
        "Minimum trace floor" => "field.trace_minimum_floor",
        "Display points" => "field.display_points",
        "Hide below mole fraction (display only)" => "field.display_cutoff",
        "Candidate temperature lower [K]" => "candidate.bounds",
        "Candidate temperature upper [K]" => "candidate.bounds",
        "Max candidates" => "candidate.max",
        "Temperature [K]" | "Start [K]" | "End [K]" | "Solved points" | "Lower bound [K]"
        | "Upper bound [K]" | "Initial seed [K]" => "temperature.range",
        _ => "field.generic",
    };
    tooltip(language, key)
}

/// Short explanations for selectable values inside segmented controls.
/// Returns an option tooltip from the structured resources.
pub fn option_tooltip(language: EquilibriumHelpLanguage, key: &str) -> &'static str {
    resource_tooltip(language, key).unwrap_or("Equilibrium selection.")
}

/// Temporary compatibility lookup for option keys not yet represented by a
/// structured resource. New UI code must use [`option_tooltip`].
#[deprecated(note = "use option_tooltip; remaining legacy branches are being removed")]
pub fn legacy_option_tooltip(language: EquilibriumHelpLanguage, key: &str) -> &'static str {
    if let Some(resource_text) = resource_tooltip(language, key) {
        return resource_text;
    }
    let english = match key {
        "state.gas" => "Gas phase; ideal-gas composition and pressure activities apply.",
        "state.liquid" => "Liquid phase candidate.",
        "state.solid" => "Solid phase candidate.",
        "model.ideal_gas" => "Ideal gas thermodynamic model.",
        "model.ideal_solution" => "Ideal solution model for a mixed phase.",
        "model.pure_condensed" => "Pure condensed-phase model.",
        "matching.subset" => {
            "Allow records containing the requested elements plus additional elements."
        }
        "matching.exact" => "Require the record element set to equal the requested set.",
        "basis.moles" => "Plot component amounts in moles.",
        "basis.fractions" => "Plot component mole fractions.",
        "basis.phase_totals" => "Plot total amount per phase.",
        "scale.linear" => "Display values on their ordinary linear scale.",
        "scale.log" => "Display positive values on a base-10 logarithmic scale.",
        "interpolation.linear" => "Interpolate display values linearly.",
        "interpolation.log" => "Interpolate positive display values in log space.",
        _ => "Equilibrium selection.",
    };
    if language == EquilibriumHelpLanguage::English {
        return english;
    }
    match key {
        "state.gas" => "Газовая фаза; применяются идеальный газ и активности давления.",
        "state.liquid" => "Кандидат жидкой фазы.",
        "state.solid" => "Кандидат твердой фазы.",
        "model.ideal_gas" => "Идеальная газовая модель.",
        "model.ideal_solution" => "Идеальный раствор для смешанной фазы.",
        "model.pure_condensed" => "Модель чистой конденсированной фазы.",
        "matching.subset" => "Разрешить дополнительные элементы сверх запрошенных.",
        "matching.exact" => "Требовать точное совпадение множества элементов.",
        "basis.moles" => "График количеств компонентов в молях.",
        "basis.fractions" => "График мольных долей компонентов.",
        "basis.phase_totals" => "График общего количества по фазам.",
        "scale.linear" => "Обычный линейный масштаб отображения.",
        "scale.log" => "Логарифмический масштаб по основанию 10 для положительных значений.",
        "interpolation.linear" => "Линейная интерполяция отображаемых значений.",
        "interpolation.log" => "Логарифмическая интерполяция положительных значений.",
        _ => english,
    }
}

pub fn tab_key(tab: &str) -> &'static str {
    match tab {
        "Setup" => "Setup",
        "Phase control" => "Phase control",
        "Libraries" => "Libraries",
        "Numerics" => "Numerics",
        "Output" => "Output",
        "Results" => "Results",
        _ => "Unknown",
    }
}

#[cfg(test)]
mod tests {
    use super::{
        EquilibriumHelpLanguage, section_from_document, tab_key, text, text_with_fallback, tooltip,
        tooltip_strict,
    };

    #[test]
    fn every_equilibrium_tab_has_bilingual_help() {
        for tab in [
            "Setup",
            "Phase control",
            "Libraries",
            "Numerics",
            "Output",
            "Results",
        ] {
            let key = tab_key(tab);
            assert!(
                !text(EquilibriumHelpLanguage::English, key)
                    .trim()
                    .is_empty()
            );
            assert!(
                !text(EquilibriumHelpLanguage::Russian, key)
                    .trim()
                    .is_empty()
            );
        }
    }

    #[test]
    fn unknown_help_key_degrades_without_panicking() {
        assert_eq!(
            text(EquilibriumHelpLanguage::English, "Unknown"),
            "Help text is unavailable for this tab."
        );
        assert_eq!(
            text_with_fallback(EquilibriumHelpLanguage::Russian, "Unknown"),
            "Help text is unavailable for this tab."
        );
    }

    #[test]
    fn incomplete_localized_document_uses_english_section_without_side_effects() {
        const ENGLISH: &str = "# Help\n\n## Setup\nEnglish setup.\n\n## Output\nEnglish output.\n";
        const INCOMPLETE_RUSSIAN: &str = "# Помощь\n\n## Setup\nРусская настройка.\n";

        let localized = section_from_document(INCOMPLETE_RUSSIAN, "Setup")
            .or_else(|| section_from_document(ENGLISH, "Setup"));
        let fallback = section_from_document(INCOMPLETE_RUSSIAN, "Output")
            .or_else(|| section_from_document(ENGLISH, "Output"));
        let missing = section_from_document(INCOMPLETE_RUSSIAN, "Unknown")
            .or_else(|| section_from_document(ENGLISH, "Unknown"));

        assert_eq!(localized, Some("Русская настройка.\n"));
        assert_eq!(fallback, Some("English output.\n"));
        assert_eq!(missing, None);
    }

    #[test]
    fn embedded_resources_contain_the_same_required_tab_sections() {
        for tab in [
            "Setup",
            "Phase control",
            "Libraries",
            "Numerics",
            "Output",
            "Results",
        ] {
            let english = text(EquilibriumHelpLanguage::English, tab);
            let russian = text(EquilibriumHelpLanguage::Russian, tab);
            assert!(
                english.lines().count() >= 2,
                "English help is too short: {tab}"
            );
            assert!(
                russian.lines().count() >= 2,
                "Russian help is too short: {tab}"
            );
        }
    }

    #[test]
    fn help_sections_contain_control_specific_guidance() {
        let required = [
            ("Setup", ["Pressure", "P,H", "phase"]),
            ("Phase control", ["hysteresis", "TPD", "iteration"]),
            ("Libraries", ["Offline", "NIST", "provenance"]),
            ("Numerics", ["backend", "tolerance", "scaling"]),
            ("Output", ["timing", "PCHIP", "display"]),
            ("Results", ["snapshot", "diagnostic", "provenance"]),
        ];
        for (tab, terms) in required {
            for language in [
                EquilibriumHelpLanguage::English,
                EquilibriumHelpLanguage::Russian,
            ] {
                let section = text(language, tab).to_lowercase();
                for term in terms {
                    assert!(
                        section.contains(&term.to_lowercase()),
                        "{language:?} help for {tab} lacks guidance for {term}"
                    );
                }
            }
        }
    }

    #[test]
    fn help_sections_define_compact_control_contracts() {
        for tab in [
            "Setup",
            "Phase control",
            "Libraries",
            "Numerics",
            "Output",
            "Results",
        ] {
            let english = text(EquilibriumHelpLanguage::English, tab);
            assert!(
                english.contains("### Control contract")
                    && english.contains("Visible when")
                    && english.contains("Change class"),
                "English help for {tab} lacks its control contract"
            );

            let russian = text(EquilibriumHelpLanguage::Russian, tab);
            assert!(
                russian.contains("### Контракт элементов управления")
                    && russian.contains("Когда виден")
                    && russian.contains("Класс изменения"),
                "Russian help for {tab} lacks its control contract"
            );
        }
    }

    #[test]
    fn help_uses_the_rendered_control_terminology() {
        let required = [
            ("Setup", ["Reference pressure [Pa]", "Model"]),
            ("Libraries", ["Engine default", "Explicit policy"]),
            ("Numerics", ["Production default cascade", "Single backend"]),
        ];
        for (tab, labels) in required {
            for language in [
                EquilibriumHelpLanguage::English,
                EquilibriumHelpLanguage::Russian,
            ] {
                let section = text(language, tab)
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ");
                for label in labels {
                    assert!(
                        section.contains(label),
                        "{language:?} help for {tab} does not use rendered label {label}"
                    );
                }
            }
        }
    }

    #[test]
    fn tooltip_catalog_has_bilingual_entries_for_primary_controls() {
        for key in [
            "tab.setup",
            "tab.phase_control",
            "tab.libraries",
            "tab.numerics",
            "tab.output",
            "tab.results",
            "help.language",
            "action.validate",
            "action.prepare",
            "action.run",
            "action.cancel",
            "action.embedded_plot",
            "action.kithe_plot",
        ] {
            assert!(!tooltip(EquilibriumHelpLanguage::English, key).is_empty());
            assert!(!tooltip(EquilibriumHelpLanguage::Russian, key).is_empty());
        }
        assert_eq!(
            tooltip(EquilibriumHelpLanguage::Russian, "unknown"),
            "Equilibrium calculator control."
        );
        for label in [
            "Pressure [Pa]",
            "Reference pressure [Pa]",
            "Target total enthalpy [J]",
            "Phase id",
            "Substance",
            "Moles",
            "Phase epsilon",
            "Creation driving force",
            "Keep driving force",
            "Tolerance override",
            "Max iterations override",
            "Display points",
            "Hide below mole fraction (display only)",
        ] {
            assert!(!super::field_tooltip(EquilibriumHelpLanguage::English, label).is_empty());
            assert!(!super::field_tooltip(EquilibriumHelpLanguage::Russian, label).is_empty());
        }
        for key in [
            "phase.initial_set",
            "solver.ph_route",
            "diagnostics.phase_trace",
            "diagnostics.range_trace",
            "output.interpolation",
            "output.pchip",
            "output.clamp",
        ] {
            assert!(!tooltip(EquilibriumHelpLanguage::English, key).is_empty());
            assert!(!tooltip(EquilibriumHelpLanguage::Russian, key).is_empty());
        }
    }

    #[test]
    fn tooltip_control_inventory_is_non_empty_in_both_languages() {
        for key in super::registered_tooltip_ids() {
            assert!(
                tooltip_strict(EquilibriumHelpLanguage::English, key).is_some(),
                "missing specific English tooltip: {key}"
            );
            assert!(
                tooltip_strict(EquilibriumHelpLanguage::Russian, key).is_some(),
                "missing specific Russian tooltip or English fallback: {key}"
            );
        }
    }

    #[test]
    fn russian_tooltip_fallbacks_are_explicit_and_registered() {
        for key in super::RUSSIAN_TOOLTIP_FALLBACK_KEYS {
            assert!(
                super::tooltip_descriptor(key).is_some(),
                "Russian fallback key is not registered: {key}"
            );
        }

        let implicit_fallbacks = super::registered_tooltip_ids()
            .filter(|key| {
                tooltip(EquilibriumHelpLanguage::English, key)
                    == tooltip(EquilibriumHelpLanguage::Russian, key)
            })
            .collect::<Vec<_>>();
        assert_eq!(
            implicit_fallbacks,
            super::RUSSIAN_TOOLTIP_FALLBACK_KEYS,
            "every untranslated Russian tooltip must be declared explicitly"
        );
    }

    #[test]
    fn strict_tooltip_lookup_rejects_unregistered_or_generic_keys() {
        assert_eq!(
            tooltip(EquilibriumHelpLanguage::English, "not-a-real-control"),
            "Equilibrium calculator control."
        );
        assert!(tooltip_strict(EquilibriumHelpLanguage::English, "not-a-real-control").is_none());
        assert!(tooltip_strict(EquilibriumHelpLanguage::Russian, "field.generic").is_none());
    }

    #[test]
    fn tooltip_control_inventory_has_unique_keys() {
        let mut keys = super::registered_tooltip_ids().collect::<Vec<_>>();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), super::EQUILIBRIUM_TOOLTIP_CATALOG.len());
    }

    #[test]
    fn descriptor_lookup_and_iterator_are_derived_from_the_catalog() {
        for descriptor in super::EQUILIBRIUM_TOOLTIP_CATALOG {
            assert_eq!(
                super::tooltip_descriptor(descriptor.key.as_str()),
                Some(descriptor)
            );
        }
        assert!(super::tooltip_descriptor("not-a-real-control").is_none());
        assert_eq!(
            super::registered_tooltip_ids().count(),
            super::EQUILIBRIUM_TOOLTIP_CATALOG.len()
        );
        let key = super::EquilibriumHelpKey::try_from("action.run").unwrap();
        assert_eq!(key.as_str(), "action.run");
        assert_eq!(
            super::tooltip_typed(EquilibriumHelpLanguage::English, key),
            tooltip(EquilibriumHelpLanguage::English, "action.run")
        );
        assert!(super::EquilibriumHelpKey::try_from("not-a-real-control").is_none());
    }

    #[test]
    fn structured_tooltip_resources_have_matching_keys_and_override_fallbacks() {
        fn keys(document: &str) -> Vec<&str> {
            document
                .lines()
                .filter_map(|line| line.strip_prefix("## "))
                .collect()
        }

        let english = keys(include_str!("../assets/equilibrium_tooltips_eng.md"));
        let russian = keys(include_str!("../assets/equilibrium_tooltips_rus.md"));
        assert!(!english.is_empty());
        let mut english_unique = english.clone();
        english_unique.sort_unstable();
        english_unique.dedup();
        assert_eq!(
            english_unique.len(),
            english.len(),
            "duplicate English structured tooltip key"
        );
        let mut russian_unique = russian.clone();
        russian_unique.sort_unstable();
        russian_unique.dedup();
        assert_eq!(
            russian_unique.len(),
            russian.len(),
            "duplicate Russian structured tooltip key"
        );
        assert_eq!(
            english_unique, russian_unique,
            "English and Russian tooltip keys diverged"
        );
        for key in &english_unique {
            assert!(
                super::resource_tooltip(EquilibriumHelpLanguage::English, key).is_some(),
                "catalog key {key} has no English structured tooltip"
            );
            assert!(
                super::resource_tooltip(EquilibriumHelpLanguage::Russian, key).is_some(),
                "catalog key {key} has no Russian structured tooltip"
            );
        }
        for key in english_unique {
            assert!(
                super::tooltip_descriptor(key).is_some(),
                "uncatalogued key: {key}"
            );
            let en = tooltip(EquilibriumHelpLanguage::English, key);
            let ru = tooltip(EquilibriumHelpLanguage::Russian, key);
            assert!(!en.is_empty());
            assert!(!ru.is_empty());
            assert_ne!(en, "Equilibrium calculator control.");
            assert_ne!(ru, en, "Russian resource did not localize {key}");
        }
        assert_eq!(
            super::option_tooltip(EquilibriumHelpLanguage::English, "state.gas"),
            "Gas phase; ideal-gas composition and pressure activities apply."
        );
        assert_ne!(
            super::option_tooltip(EquilibriumHelpLanguage::Russian, "state.gas"),
            super::option_tooltip(EquilibriumHelpLanguage::English, "state.gas")
        );
    }

    #[test]
    fn catalog_registers_actions_and_field_descriptions() {
        assert!(super::registered_tooltip_ids().any(|key| key == "action.run"));
        assert!(
            super::registered_tooltip_ids().any(|key| key == "field.pressure"),
            "field-level descriptions remain catalog entries even when their widget\n             coverage belongs exclusively to story tests"
        );
    }
}
