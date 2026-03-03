use crate::Kinetics::experimental_kinetics::exp_engine_api::Ranges;

use crate::Kinetics::experimental_kinetics::experiment_series_main::TGASeries;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADomainError;

use crate::gui::experimental_kinetics_gui::interaction::InteractionState;
use crate::gui::experimental_kinetics_gui::settings::Settings;
use log::info;
use simplelog::{Config, LevelFilter, SimpleLogger};

use std::sync::Once;
use tabled::{Table, Tabled};
/// Константа, определяющая диапазон отображения по умолчанию (от 0 до 10)
const DEFAULT_VIEW_RANGE: (f64, f64) = (0.0, 10.0);

//======================================================================
// ERRORS
#[derive(Debug)]
pub enum TGAGUIError {
    TGADomainError(TGADomainError),
    SettingsErrors(String),
    BindingError(String),
}

impl From<TGADomainError> for TGAGUIError {
    fn from(err: TGADomainError) -> Self {
        TGAGUIError::TGADomainError(err)
    }
}

fn parse_log_level(level: Option<&str>) -> LevelFilter {
    match level.unwrap_or("info").trim().to_ascii_lowercase().as_str() {
        "off" => LevelFilter::Off,
        "error" => LevelFilter::Error,
        "warn" | "warning" => LevelFilter::Warn,
        "info" => LevelFilter::Info,
        "debug" => LevelFilter::Debug,
        "trace" => LevelFilter::Trace,
        _ => LevelFilter::Info,
    }
}

fn init_logging(settings: &Settings) {
    static LOGGER_INIT: Once = Once::new();
    let level = parse_log_level(settings.log_level());
    let requested_level = settings.log_level().unwrap_or("info").to_string();

    LOGGER_INIT.call_once(|| {
        let _ = SimpleLogger::init(level, Config::default());
    });

    log::set_max_level(level);
    info!(
        "Experimental kinetics logger initialized with level '{}' ({:?})",
        requested_level, level
    );
}
//======================================================================
/// TGA Plot Curve - represents experimental data curve
#[derive(Debug, Clone)]
pub struct PlotCurve {
    pub experiment_index: usize,
    pub experiment_id: String,
    pub plot_short_name: String,
    pub color: [u8; 3],
    pub selected: bool,
    pub points: Vec<[f64; 2]>,
    pub x_name: String,
    pub y_name: String,
    pub ranges: Ranges,
    pub highlighted: bool,
    pub shown: bool,
}

/// Типовые цвета для кривых графика / Typical preset colours for plot curves
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Colours {
    Blue,
    Red,
    Green,
    Yellow,
    Orange,
    Purple,
    Cyan,
    Magenta,
    Black,
    White,
    Gray,
}

impl Colours {
    /// Преобразует вариант enum в RGB-массив / Converts enum variant into RGB array
    pub fn as_rgb(self) -> [u8; 3] {
        match self {
            Colours::Blue => [0, 0, 255],
            Colours::Red => [255, 0, 0],
            Colours::Green => [0, 255, 0],
            Colours::Yellow => [255, 255, 0],
            Colours::Orange => [255, 165, 0],
            Colours::Purple => [128, 0, 128],
            Colours::Cyan => [0, 255, 255],
            Colours::Magenta => [255, 0, 255],
            Colours::Black => [0, 0, 0],
            Colours::White => [255, 255, 255],
            Colours::Gray => [128, 128, 128],
        }
    }
}

impl PlotCurve {
    pub fn find_nearest_distance(&self, point: [f64; 2]) -> Option<f64> {
        if !self.shown || self.points.is_empty() {
            return None;
        }

        // Normalize distances to keep picking stable for very different axis scales.
        let x_span = (self.ranges.x_max - self.ranges.x_min)
            .abs()
            .max(f64::EPSILON);
        let y_span = (self.ranges.y_max - self.ranges.y_min)
            .abs()
            .max(f64::EPSILON);

        let point_to_segment_normalized = |a: [f64; 2], b: [f64; 2], p: [f64; 2]| -> f64 {
            let abx = b[0] - a[0];
            let aby = b[1] - a[1];
            let apx = p[0] - a[0];
            let apy = p[1] - a[1];
            let ab_len2 = abx * abx + aby * aby;

            let t = if ab_len2 > 0.0 {
                (apx * abx + apy * aby) / ab_len2
            } else {
                0.0
            }
            .clamp(0.0, 1.0);

            let cx = a[0] + t * abx;
            let cy = a[1] + t * aby;

            let dx = (p[0] - cx) / x_span;
            let dy = (p[1] - cy) / y_span;
            (dx * dx + dy * dy).sqrt()
        };

        if self.points.len() == 1 {
            let dx = (point[0] - self.points[0][0]) / x_span;
            let dy = (point[1] - self.points[0][1]) / y_span;
            return Some((dx * dx + dy * dy).sqrt());
        }

        self.points
            .windows(2)
            .map(|w| point_to_segment_normalized(w[0], w[1], point))
            .min_by(|a, b| a.total_cmp(b))
    }

    pub fn set_shown(&mut self, shown: bool) {
        self.shown = shown;
    }

    pub fn is_shown(&self) -> bool {
        self.shown
    }

    pub fn get_name(&self) -> String {
        self.experiment_id.clone()
    }

    pub fn get_label(&self) -> String {
        format!("{}: {} vs {}", self.experiment_id, self.y_name, self.x_name)
    }
    pub fn get_short_name(&self) -> String {
        self.plot_short_name.clone()
    }

    /// Устанавливает цвет по enum-варианту / Sets curve colour using enum variant
    pub fn set_colour(&mut self, colour: Colours) {
        self.color = colour.as_rgb();
    }

    /// Устанавливает цвет напрямую по RGB / Sets curve colour directly from RGB
    pub fn set_colour_rgb(&mut self, rgb: [u8; 3]) {
        self.color = rgb;
    }
}

impl Default for PlotCurve {
    fn default() -> Self {
        Self {
            experiment_index: 0,
            experiment_id: String::new(),
            plot_short_name: String::new(),
            color: [0, 0, 255],

            selected: false,
            points: Vec::new(),
            x_name: String::new(),
            y_name: String::new(),
            ranges: Ranges::default(),
            highlighted: false,
            shown: true,
        }
    }
}
/// Основная модель данных приложения
///
/// Хранит все данные графиков, состояние выделения,
/// параметры отображения и другие данные приложения
#[derive(Debug, Clone)]
pub struct PlotModel {
    /// Доменные данные
    pub series: TGASeries,
    /// Вектор графиков для отображения
    pub plots: Vec<PlotCurve>,

    pub interaction: InteractionState,
    /// Количество точек для отображения каждого графика
    pub settings: Settings,
    /// Флаг, указывающий, что был запрошен сброс вида
    pub reset_view_requested: bool,
    pub plot_recreation_required: bool,
    pub message: String,
}

impl Default for PlotModel {
    /// Реализация значения по умолчанию для модели
    ///

    /// стандартными параметрами отображения и сброшенным состоянием
    fn default() -> Self {
        let settings = Settings::new().unwrap_or_else(|_| Settings::default());
        init_logging(&settings);

        // Print settings as a table to terminal
        print_settings_table(&settings);

        Self {
            series: TGASeries::new(),
            plots: Vec::new(),

            interaction: InteractionState::default(),
            settings,
            reset_view_requested: false,
            plot_recreation_required: false,
            message: String::new(),
        }
    }
}

/// Helper struct for displaying settings in a table format
#[derive(Tabled)]
struct SettingsRow {
    #[tabled(rename = "Parameter")]
    parameter: String,
    #[tabled(rename = "Value")]
    value: String,
}

/// Prints settings as a formatted table to the terminal
fn print_settings_table(settings: &Settings) {
    let mut rows = Vec::new();

    // Add calibration line settings
    if let Some(cal) = settings.calibration_line() {
        rows.push(SettingsRow {
            parameter: "Calibration k".to_string(),
            value: format!("{:.4}", cal.k()),
        });
        rows.push(SettingsRow {
            parameter: "Calibration b".to_string(),
            value: format!("{:.4}", cal.b()),
        });
    }

    // Add number of points
    if let Some(n) = settings.n_points() {
        rows.push(SettingsRow {
            parameter: "Number of points".to_string(),
            value: n.to_string(),
        });
    }

    // Add symbolic expression
    if let Some(expr) = settings.symbolic_expression() {
        rows.push(SettingsRow {
            parameter: "Symbolic expression".to_string(),
            value: expr.to_string(),
        });
    }

    // Add log level
    if let Some(log) = settings.log_level() {
        rows.push(SettingsRow {
            parameter: "Log level".to_string(),
            value: log.to_string(),
        });
    }

    // Create and print the table
    let table = Table::new(rows).to_string();
    println!("\n=== TGA Application Settings ===");
    println!("{}", table);
    println!("================================\n");
}

//============================================================================================
impl PlotModel {
    /// Создает новую модель с параметрами по умолчанию
    ///
    /// # Возвращает
    /// Новая модель данных
    pub fn new() -> Self {
        Self::default()
    }

    /// Выбирает график по индексу
    ///
    /// Снимает выбор со всех графиков и устанавливает выбор
    /// только на график с указанным индексом
    ///
    /// # Параметры
    /// * `index` - индекс графика для выбора
    pub fn select_curve(&mut self, index: usize) {
        // Проходим по всем графикам и устанавливаем выбор только на нужный
        for (i, curve) in self.plots.iter_mut().enumerate() {
            curve.selected = i == index;
        }
    }
    /// Снимает выбор со всех графиков
    pub fn clear_selection(&mut self) {
        // Проходим по всем графикам и снимаем выбор
        for curve in &mut self.plots {
            curve.selected = false;
        }
    }
    /// Возвращает индекс выбранного графика
    ///
    /// # Возвращает
    /// Some(индекс) если есть выбранный график, None если нет
    pub fn get_selected_curve_index(&self) -> Option<usize> {
        self.plots.iter().position(|c| c.selected)
    }

    pub fn get_experiment_by_selected_curve(&self) -> Result<String, TGAGUIError> {
        if let Some(number) = self.get_selected_curve_index() {
            return Ok(self.plots[number].experiment_id.clone());
        } else if self.is_only_one_shown() {
            let id = self
                .plots
                .iter()
                .find(|curve| curve.is_shown())
                .map(|curve| curve.experiment_id.clone())
                .unwrap_or_default();
            if !id.is_empty() {
                return Ok(id);
            }
        }

        Err(TGAGUIError::BindingError(
            "No selected curve. Select a curve first.".to_string(),
        ))
    }

    /// Returns true when all currently shown curves belong to one experiment.
    pub fn is_only_one_shown(&self) -> bool {
        let mut shown_curves = self.plots.iter().filter(|curve| curve.is_shown());
        let Some(first) = shown_curves.next() else {
            return false;
        };
        shown_curves.all(|curve| curve.experiment_id == first.experiment_id)
    }
    //======================================================================================
    // INTERACTION THIN WRAPPER
    //=====================================================================================
    /// Начинает процесс выделения области
    ///
    /// # Параметры
    /// * `start_point` - начальная точка выделения [x, y]
    pub fn start_selection(&mut self, start_point: [f64; 2]) {
        self.interaction.start_selection(start_point);
        // Очищаем предыдущие подсветки
        self.clear_highlights();
    }

    /// Обновляет процесс панорамирования
    ///
    /// # Параметры
    /// * `current_point` - текущая позиция панорамирования [x, y]
    pub fn update_pan(&mut self, current_point: [f64; 2]) {
        self.interaction.update_pan(current_point);
    }

    /// Обновляет процесс выделения
    ///
    /// # Параметры
    /// * `end_point` - конечная точка выделения [x, y]
    pub fn update_selection(&mut self, end_point: [f64; 2]) {
        if self.interaction.update_selection(end_point) {
            self.update_highlighted_curves();
        }
    }
    /// Завершает процесс выделения
    pub fn end_selection(&mut self) {
        self.interaction.end_selection();
        // Обновляем подсветки в последний раз
        self.update_highlighted_curves();
    }
    /// Очищает прямоугольник выделения
    pub fn clear_selection_rect(&mut self) {
        self.interaction.clear_selection_rect();
        self.clear_highlights();
        // Запрос обновления UI
    }
    /// Масштабирует отображение
    ///
    /// # Параметры
    /// * `factor` - коэффициент масштабирования
    /// * `center` - центральная точка масштабирования [x, y]
    pub fn zoom(&mut self, factor: f64, center: [f64; 2]) {
        self.interaction.zoom(factor, center);
    }
    /// Масштабирует отображение под выделенную область
    pub fn fit_to_selection(&mut self) {
        self.interaction.fit_to_selection();
    }
    //==============================================================================
    /// Обновляет статус подсветки точек графиков
    ///
    /// Проверяет, какие точки графиков находятся внутри
    /// текущего прямоугольника выделения и устанавливает
    /// соответствующий статус подсветки
    pub fn update_highlighted_curves(&mut self) {
        if let Some(rect) = &self.interaction.selection_rect {
            let (x_min, x_max, y_min, y_max) = rect.bounds();
            for curve in &mut self.plots {
                let has_point_in_selection = curve
                    .points
                    .iter()
                    .any(|&[x, y]| x >= x_min && x <= x_max && y >= y_min && y <= y_max);
                // Устанавливаем статус подсветки для графика
                curve.highlighted = has_point_in_selection;
            }
        } else {
            // Очищаем подсветки, если нет прямоугольника выделения
            self.clear_highlights();
        }
    }
    /// Очищает подсветку всех графиков
    pub fn clear_highlights(&mut self) {
        // Проходим по всем графикам и снимаем подсветку
        for curve in &mut self.plots {
            curve.highlighted = false;
        }
    }

    /// Сбрасывает отображение к значениям по умолчанию
    /// Now calculates ranges from all visible plots instead of using hardcoded values
    pub fn reset_view(&mut self) {
        // Calculate ranges from all visible plots
        let visible_plots: Vec<_> = self.plots.iter().filter(|p| p.is_shown()).collect();

        if visible_plots.is_empty() {
            // No visible plots - use default ranges
            self.interaction.view_range = DEFAULT_VIEW_RANGE;
            self.interaction.view_y_range = (-1.0, 1.0);
        } else {
            // Calculate combined ranges from all visible plots
            let mut x_min = f64::INFINITY;
            let mut x_max = f64::NEG_INFINITY;
            let mut y_min = f64::INFINITY;
            let mut y_max = f64::NEG_INFINITY;

            for plot in visible_plots {
                x_min = x_min.min(plot.ranges.x_min);
                x_max = x_max.max(plot.ranges.x_max);
                y_min = y_min.min(plot.ranges.y_min);
                y_max = y_max.max(plot.ranges.y_max);
            }

            // Add margin and ensure non-degenerate bounds even for flat/single-point curves.
            let x_span = (x_max - x_min).abs();
            let y_span = (y_max - y_min).abs();
            let x_floor = x_max.abs().max(x_min.abs()) * 1e-6 + 1e-9;
            let y_floor = y_max.abs().max(y_min.abs()) * 1e-6 + 1e-9;
            let x_margin = (x_span * 0.05).max(x_floor);
            let y_margin = (y_span * 0.05).max(y_floor);

            self.interaction.view_range = (x_min - x_margin, x_max + x_margin);
            self.interaction.view_y_range = (y_min - y_margin, y_max + y_margin);
        }

        self.reset_view_requested = true;
    }
    /// Получает и сбрасывает флаг запроса сброса вида
    ///
    /// # Возвращает
    /// Предыдущее значение флага запроса сброса
    pub fn take_reset_view_request(&mut self) -> bool {
        let was = self.reset_view_requested;
        self.reset_view_requested = false;
        was
    }
    /// Находит ближайший график к заданной точке
    ///
    /// # Параметры
    /// * `point` - точка [x, y] для поиска ближайшего графика
    ///
    /// # Возвращает
    /// Some(индекс) если найден близкий график, None если нет
    pub fn find_nearest_curve(&self, point: [f64; 2]) -> Option<usize> {
        // Normalized (dimensionless) pick radius.
        let tolerance = 0.06;
        self.plots
            .iter()
            .enumerate()
            .filter(|(_, curve)| curve.is_shown())
            .filter_map(|(i, curve)| {
                curve
                    .find_nearest_distance(point)
                    .filter(|&d| d < tolerance)
                    .map(|d| (i, d))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _)| i)
    }

    /*
     pub fn reset_view(&mut self) {
        self.interaction.view_range = DEFAULT_VIEW_RANGE; // Сбрасываем диапазон по X
        // Сбрасываем диапазон по Y на основе текущих амплитуд
        let y_extent = (self.max_amplitude() * 1.5).max(1.0);
        self.interaction.view_y_range = (-y_extent, y_extent);
        self.reset_view_requested = true; // Устанавливаем флаг запроса сброса
    }
     */
}
