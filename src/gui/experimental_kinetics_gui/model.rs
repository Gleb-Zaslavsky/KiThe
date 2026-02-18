use crate::Kinetics::experimental_kinetics::exp_engine_api::{PlotSeries, Ranges, ViewRange};
use crate::Kinetics::experimental_kinetics::experiment_series_main::{TGAExperiment, TGASeries};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADomainError, Unit};
use crate::gui::experimental_kinetics_gui::interaction::InteractionState;
use crate::gui::experimental_kinetics_gui::settings::Settings;
use std::path::Path;
/// Константа, определяющая диапазон отображения по умолчанию (от 0 до 10)
const DEFAULT_VIEW_RANGE: (f64, f64) = (0.0, 10.0);

//======================================================================
// ERRORS
#[derive(Debug)]
pub enum TGAGUIError {
    TGADomainError(TGADomainError),
    SettingsErrors(String),
}

impl From<TGADomainError> for TGAGUIError {
    fn from(err: TGADomainError) -> Self {
        TGAGUIError::TGADomainError(err)
    }
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
    pub n_points: usize,
    /// Флаг, указывающий, что был запрошен сброс вида
    pub reset_view_requested: bool,
}

impl Default for PlotModel {
    /// Реализация значения по умолчанию для модели
    ///
    /// Создает модель с двумя графиками (синус и косинус),
    /// стандартными параметрами отображения и сброшенным состоянием
    fn default() -> Self {
        let settings = Settings::new().unwrap_or_else(|_| Settings::default());
        let n_points = settings.n_points().unwrap_or(1000);

        Self {
            series: TGASeries::new(),
            plots: Vec::new(),

            interaction: InteractionState::default(),
            n_points,
            reset_view_requested: false,
        }
    }
}
//===========================================================================================
impl PlotModel {
    //=================================================================
    // SETTERS
    //=======================================================================
    pub fn push_from_file(&mut self, path: &Path) -> Result<(), TGAGUIError> {
        self.series.push_from_file(path)?;
        Ok(())
    }
    /// list of experiments
    pub fn list_of_experiments(&self) -> Vec<String> {
        self.series.ids()
    }

    pub fn get_list_of_complex_id(&self) -> Vec<String> {
        self.series.get_list_of_complex_id()
    }

    pub fn get_experiment_by_id_mut(
        &mut self,
        id: &str,
    ) -> Result<&mut TGAExperiment, TGAGUIError> {
        let exp = self.series.get_experiment_by_id_mut(id)?;
        Ok(exp)
    }

    pub fn bind_time_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_time(id, col, unit)?;
        Ok(())
    }

    pub fn bind_temperature_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_temperature(id, col, unit)?;
        Ok(())
    }

    pub fn bind_mass_of_experiment(
        &mut self,
        id: &str,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        self.series.bind_mass(id, col, unit)?;
        Ok(())
    }

    pub fn bind_time_of_selected(&mut self, col: &str, unit: Unit) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.bind_time_of_experiment(&id, col, unit)?;
        Ok(())
    }

    pub fn bind_temperature_of_selected(
        &mut self,
        col: &str,
        unit: Unit,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.bind_temperature_of_experiment(&id, col, unit)?;
        Ok(())
    }

    pub fn bind_mass_of_selected(&mut self, col: &str, unit: Unit) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.bind_mass_of_experiment(&id, col, unit)?;
        Ok(())
    }
    //========================================================================

    //========================================================================
    // GENERIC METHODS
    pub fn apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&TGAExperiment) -> R,
    {
        let res = self.series.apply_by_id(id, op)?;
        Ok(res)
    }

    /// try_apply_by_id API method in the experiment/series facade.
    /// Calls a fallible read-only operation on one experiment selected by id.
    /// generic argument:  FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_apply_by_id<R, F>(&self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&TGAExperiment) -> Result<R, TGADomainError>,
    {
        let res = self.series.try_apply_by_id(id, op)?;
        Ok(res)
    }

    /// mutate_by_id API method in the experiment/series facade.
    /// Calls an infallible mutable operation on one experiment selected by id.
    ///  generic argument:  FnOnce(&mut TGAExperiment) -> R,
    pub fn mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&mut TGAExperiment) -> R,
    {
        let res = self.series.mutate_by_id(id, op)?;
        Ok(res)
    }

    /// try_mutate_by_id API method in the experiment/series facade.
    /// Calls a fallible mutable operation on one experiment selected by id.
    ///generic argument:  FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    pub fn try_mutate_by_id<R, F>(&mut self, id: &str, op: F) -> Result<R, TGAGUIError>
    where
        F: FnOnce(&mut TGAExperiment) -> Result<R, TGADomainError>,
    {
        let res = self.series.try_mutate_by_id(id, op)?;
        Ok(res)
    }

    /// transform_by_id API method in the experiment/series facade.
    /// Applies an infallible consuming transformation (`TGAExperiment -> TGAExperiment`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> TGAExperiment,
    pub fn transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGAGUIError>
    where
        F: FnOnce(TGAExperiment) -> TGAExperiment,
    {
        self.series.transform_by_id(id, op)?;
        Ok(())
    }

    /// try_transform_by_id API method in the experiment/series facade.
    /// Applies a fallible consuming transformation
    /// (`TGAExperiment -> Result<TGAExperiment, TGADomainError>`)
    /// to one experiment selected by id.
    /// generic argument: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    pub fn try_transform_by_id<F>(&mut self, id: &str, op: F) -> Result<(), TGAGUIError>
    where
        F: FnOnce(TGAExperiment) -> Result<TGAExperiment, TGADomainError>,
    {
        self.series.try_transform_by_id(id, op)?;

        Ok(())
    }
    //=================================================================================
    // CREATING PLOT
    //=========================================================================
    ///
    pub fn set_x(&mut self, id: &str, x_col: &str) -> Result<(), TGAGUIError> {
        self.series.set_oneframeplot_x(id, x_col)?;
        Ok(())
    }
    pub fn set_y(&mut self, id: &str, y_col: &str) -> Result<(), TGAGUIError> {
        self.series.set_oneframeplot_y(id, y_col)?;
        Ok(())
    }
    pub fn list_of_columns(&self, id: &str) -> Result<Vec<String>, TGAGUIError> {
        let list = self.series.list_of_columns(id)?;
        Ok(list)
    }
    /// max and min of x and y
    pub fn plot_range(&mut self, id: &str) -> Result<Ranges, TGAGUIError> {
        println!("\n setting plot named {} the  range!", id);
        let ranges = self.series.plot_xy_ranges(id)?;
        Ok(ranges)
    }
    // initial ranges
    pub fn give_ranges_to_plot(&mut self, id: &str) -> Result<Ranges, TGAGUIError> {
        println!("\n giving ranges to plot named {}", id);

        let ranges = self.plot_range(id)?;
        let x_min = ranges.x_min;
        let y_min = ranges.y_min;
        let x_max = ranges.x_max;
        let y_max = ranges.y_max;
        self.interaction.view_range = (x_min, x_max);
        self.interaction.view_y_range = (y_min, y_max);
        Ok(ranges)
    }
    pub fn calculate_samplepoints_for_plot(
        &mut self,
        id: &str,
        x_range: [f64; 2],
    ) -> Result<PlotSeries, TGAGUIError> {
        let range = Some({
            ViewRange {
                t_min: x_range[0],
                t_max: x_range[1],
            }
        });
        let plot = self.series.sample_oneframeplot(id, range, self.n_points)?;
        Ok(plot)
    }
    /// get number of plot by itd String id
    pub fn index_by_id(&self, id: &str) -> Result<usize, TGAGUIError> {
        let res = self.series.index_by_id(id)?;
        Ok(res)
    }
    pub fn create_points_for_curve(&mut self, id: &str) -> Result<(), TGAGUIError> {
        println!("calculating points for curve");
        let ranges = self.give_ranges_to_plot(id)?;
        let x0 = self.interaction.view_range.0;
        let x1 = self.interaction.view_range.1;
        let x_range = [x0, x1];
        let plot = self.calculate_samplepoints_for_plot(id, x_range)?;
        let index = self.index_by_id(id)?;
        self.push_to_vec_of_curves(plot, id, index, ranges)?;
        Ok(())
    }

    pub fn create_plot_short_name(&self, id: &str, x_name: &String, y_name: &String) -> String {
        format!("{}_{}_{}", id.to_string(), x_name, y_name)
    }

    pub fn push_to_vec_of_curves(
        &mut self,
        plot: PlotSeries,
        id: &str,
        index: usize,
        ranges: Ranges,
    ) -> Result<(), TGAGUIError> {
        let mut new_curve = PlotCurve::default();
        new_curve.experiment_id = id.to_string().clone();
        new_curve.experiment_index = index;
        let x_name = plot.name_x;
        let y_name = plot.name_y;
        let short_name = self.create_plot_short_name(id, &x_name, &y_name);
        new_curve.x_name = x_name;
        new_curve.y_name = y_name;
        let x_vec = plot.x;
        let y_vec = plot.y;
        let mut points = Vec::new();
        for (xi, yi) in x_vec.iter().zip(y_vec) {
            points.push([*xi, yi]);
        }
        new_curve.ranges = ranges;
        new_curve.points = points;
        new_curve.plot_short_name = short_name;
        self.plots.push(new_curve);
        Ok(())
    }
    //=================================================================================
    //COLUMNS TRASFORMATION
    //================================================================================

    // column transformations
    pub fn move_time_to_zero_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.move_time_to_zero(id)?;
        Ok(())
    }

    pub fn move_time_to_zero_of_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.move_time_to_zero_for_experiment(&id)?;
        Ok(())
    }

    pub fn from_C_to_K_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.celsius_to_kelvin(id)?;
        Ok(())
    }

    pub fn from_C_to_K_of_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.from_C_to_K_for_experiment(&id)?;
        Ok(())
    }

    pub fn from_s_to_h_for_experiment(&mut self, id: &str) -> Result<(), TGAGUIError> {
        self.series.seconds_to_hours(id)?;
        Ok(())
    }

    pub fn from_s_to_h_of_selected(&mut self) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.from_s_to_h_for_experiment(&id)?;
        Ok(())
    }

    pub fn offset_column(
        &mut self,
        id: &str,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.series.offset_column(id, column, offset)?;
        Ok(())
    }
    pub fn offset_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.series.offset_column(&id, column, offset)?;
        Ok(())
    }
    pub fn add_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_column_for_selected(column, offset)
    }
    pub fn sub_column_for_selected(
        &mut self,
        column: &str,
        offset: f64,
    ) -> Result<(), TGAGUIError> {
        self.offset_column_for_selected(column, -offset)
    }

    pub fn scale_columns(
        &mut self,
        id: &str,
        cols: &[&str],
        factor: f64,
    ) -> Result<(), TGAGUIError> {
        self.series.scale_columns(id, cols, factor)?;
        Ok(())
    }

    pub fn scale_column(&mut self, id: &str, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        self.scale_columns(id, &[col], factor)?;
        Ok(())
    }

    pub fn scale_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.scale_column(&id, col, factor)?;
        Ok(())
    }

    pub fn mul_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        self.scale_column_of_selected(col, factor)
    }

    pub fn div_column_of_selected(&mut self, col: &str, factor: f64) -> Result<(), TGAGUIError> {
        self.scale_column_of_selected(col, 1.0 / factor)
    }

    pub fn exp_column(&mut self, id: &str, col_name: &str) -> Result<(), TGAGUIError> {
        self.series.exp_column(id, col_name)?;
        Ok(())
    }

    pub fn exp_column_for_selected(&mut self, col_name: &str) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.exp_column(&id, col_name)?;
        Ok(())
    }

    pub fn ln_column(&mut self, id: &str, col_name: &str) -> Result<(), TGAGUIError> {
        self.series.ln_column(id, col_name)?;
        Ok(())
    }

    pub fn ln_column_for_selected(&mut self, col_name: &str) -> Result<(), TGAGUIError> {
        let id = self.get_experiment_by_selected_curve()?;
        self.ln_column(&id, col_name)?;
        Ok(())
    }
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
        let id = if let Some(number) = self.get_selected_curve_index() {
            self.plots[number].experiment_id.clone()
        } else if self.is_only_one_shown() {
            self.plots
                .iter()
                .find(|curve| curve.is_shown())
                .map(|curve| curve.experiment_id.clone())
                .unwrap_or_default()
        } else {
            String::new()
        };
        Ok(id)
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
    /// * `current_x` - текущая x-координата панорамирования
    pub fn update_pan(&mut self, current_x: f64) {
        self.interaction.update_pan(current_x);
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
    /// * `center` - центральная точка масштабирования
    pub fn zoom(&mut self, factor: f64, center: f64) {
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
    pub fn reset_view(&mut self) {
        self.interaction.view_range = DEFAULT_VIEW_RANGE;
        self.interaction.view_y_range = (-1.0, 1.0);
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
        let tolerance = 0.03;
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
