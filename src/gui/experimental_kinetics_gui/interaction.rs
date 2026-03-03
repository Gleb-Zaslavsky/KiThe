/// Константа, определяющая диапазон отображения по умолчанию (от 0 до 10)
const DEFAULT_VIEW_RANGE: (f64, f64) = (0.0, 10.0);

/// Структура прямоугольника выделения
///
/// Хранит координаты начальной и конечной точек прямоугольника выделения
#[derive(Debug, Clone, Copy)]
pub struct SelectionRect {
    /// Начальная точка выделения [x, y]
    pub start: [f64; 2],
    /// Конечная точка выделения [x, y]
    pub end: [f64; 2],
}

impl SelectionRect {
    /// Создает новый прямоугольник выделения
    ///
    /// # Параметры
    /// * `start` - начальная точка [x, y]
    /// * `end` - конечная точка [x, y]
    ///
    /// # Возвращает
    /// Новый экземпляр SelectionRect
    pub fn new(start: [f64; 2], end: [f64; 2]) -> Self {
        Self { start, end }
    }

    /// Проверяет, активен ли прямоугольник выделения
    ///
    /// # Возвращает
    /// true, если прямоугольник активен (начальная и конечная точки различны)
    pub fn is_active(&self) -> bool {
        self.start != self.end
    }

    /// Нормализует прямоугольник выделения
    ///
    /// Упорядочивает координаты так, чтобы начальная точка
    /// была в левом нижнем углу, а конечная - в правом верхнем
    ///
    /// # Возвращает
    /// Новый нормализованный прямоугольник выделения
    pub fn normalized(&self) -> Self {
        // Находим минимальные и максимальные координаты
        let min_x = self.start[0].min(self.end[0]);
        let max_x = self.start[0].max(self.end[0]);
        let min_y = self.start[1].min(self.end[1]);
        let max_y = self.start[1].max(self.end[1]);

        Self {
            start: [min_x, min_y],
            end: [max_x, max_y],
        }
    }

    /// Проверяет, содержит ли прямоугольник заданную точку
    ///
    /// # Параметры
    /// * `point` - точка [x, y] для проверки
    ///
    /// # Возвращает
    /// true, если точка находится внутри прямоугольника
    #[allow(dead_code)]
    pub fn contains_point(&self, point: [f64; 2]) -> bool {
        // Нормализуем прямоугольник для корректной проверки
        let normalized = self.normalized();
        // Проверяем, находится ли точка внутри прямоугольника
        point[0] >= normalized.start[0]
            && point[0] <= normalized.end[0]
            && point[1] >= normalized.start[1]
            && point[1] <= normalized.end[1]
    }

    /// Возвращает границы прямоугольника
    ///
    /// # Возвращает
    /// Кортеж (x_min, x_max, y_min, y_max) с границами прямоугольника
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        // Нормализуем прямоугольник для получения корректных границ
        let normalized = self.normalized();
        (
            normalized.start[0], // x_min
            normalized.end[0],   // x_max
            normalized.start[1], // y_min
            normalized.end[1],   // y_max
        )
    }
}
#[derive(Debug, Clone)]
pub struct InteractionState {
    /// Прямоугольник текущего выделения (если есть)
    pub selection_rect: Option<SelectionRect>,
    /// Флаг, указывающий, выполняется ли выделение
    pub is_selecting: bool,
    /// Начальная точка выделения (если выполняется выделение)
    pub selection_start: Option<[f64; 2]>,
    /// Флаг, указывающий, выполняется ли панорамирование
    pub is_panning: bool,
    /// Последняя позиция панорамирования (для вычисления смещения)
    pub last_pan_pos: Option<[f64; 2]>,
    /// Диапазон отображения по оси X (минимум и максимум)
    pub view_range: (f64, f64),
    /// Диапазон отображения по оси Y (минимум и максимум)
    pub view_y_range: (f64, f64),
}

impl InteractionState {
    pub fn default() -> Self {
        Self {
            selection_rect: None,           // Нет активного выделения
            is_selecting: false,            // Выделение не выполняется
            selection_start: None,          // Нет начальной точки выделения
            is_panning: false,              // Панорамирование не выполняется
            last_pan_pos: None,             // Нет последней позиции панорамирования
            view_range: DEFAULT_VIEW_RANGE, // Диапазон отображения по умолчанию
            view_y_range: (-1.5, 1.5),      // Диапазон отображения по Y
        }
    }
    /// Сдвигает отображение по оси X
    ///
    /// # Параметры
    /// * `delta` - величина сдвига
    pub fn pan(&mut self, delta: f64) {
        self.view_range.0 += delta; // Сдвигаем минимальную границу
        self.view_range.1 += delta; // Сдвигаем максимальную границу
    }

    pub fn pan_xy(&mut self, delta_x: f64, delta_y: f64) {
        self.view_range.0 += delta_x;
        self.view_range.1 += delta_x;
        self.view_y_range.0 += delta_y;
        self.view_y_range.1 += delta_y;
    }
    /// Начинает процесс выделения области
    ///
    /// # Параметры
    /// * `start_point` - начальная точка выделения [x, y]
    pub fn start_selection(&mut self, start_point: [f64; 2]) {
        self.is_selecting = true; // Устанавливаем флаг выделения
        self.selection_start = Some(start_point); // Сохраняем начальную точку
        // Создаем прямоугольник выделения с начальной точкой как началом и концом
        self.selection_rect = Some(SelectionRect::new(start_point, start_point));
    }

    /// Начинает процесс панорамирования
    ///
    /// # Параметры
    /// * `start_x` - начальная x-координата панорамирования
    pub fn start_pan(&mut self, start_point: [f64; 2]) {
        self.is_panning = true; // Устанавливаем флаг панорамирования
        self.last_pan_pos = Some(start_point); // Сохраняем начальную позицию
    }

    /// Обновляет процесс панорамирования
    ///
    /// # Параметры
    /// * `current_x` - текущая x-координата панорамирования
    pub fn update_pan(&mut self, current_point: [f64; 2]) {
        // Проверяем, есть ли предыдущая позиция
        if let Some(prev) = self.last_pan_pos {
            // Вычисляем смещение
            let dx = current_point[0] - prev[0];
            let dy = current_point[1] - prev[1];
            // Перемещаем вид в противоположном направлении движения указателя
            self.pan_xy(-dx, -dy);
            // Обновляем последнюю позицию
            self.last_pan_pos = Some(current_point);
        }
    }

    /// Завершает процесс панорамирования
    pub fn end_pan(&mut self) {
        self.is_panning = false; // Снимаем флаг панорамирования
        self.last_pan_pos = None; // Очищаем последнюю позицию
    }

    /// Обновляет процесс выделения
    ///
    /// # Параметры
    /// * `end_point` - конечная точка выделения [x, y]
    pub fn update_selection(&mut self, end_point: [f64; 2]) -> bool {
        // Проверяем, есть ли начальная точка выделения
        if let Some(start) = self.selection_start {
            // Обновляем прямоугольник выделения
            self.selection_rect = Some(SelectionRect::new(start, end_point));
            true
        } else {
            false
        }
    }

    /// Завершает процесс выделения
    pub fn end_selection(&mut self) {
        self.is_selecting = false; // Снимаем флаг выделения
    }
    /// Очищает прямоугольник выделения
    pub fn clear_selection_rect(&mut self) {
        self.selection_rect = None; // Очищаем прямоугольник выделения
    }

    /// Масштабирует отображение
    ///
    /// # Параметры
    /// * `factor` - коэффициент масштабирования
    /// * `center` - центральная точка масштабирования
    pub fn zoom(&mut self, factor: f64, center: [f64; 2]) {
        let factor = if factor <= 0.0 { 1.0 } else { factor };
        let min_range = 1e-6_f64;

        let zoom_axis = |axis_range: (f64, f64), axis_center: f64| -> Option<(f64, f64)> {
            let mut range = axis_range.1 - axis_range.0;
            if range <= min_range {
                range = min_range;
            }

            let new_range = (range / factor).max(min_range);
            let center_used = if axis_center.is_finite()
                && axis_center >= axis_range.0
                && axis_center <= axis_range.1
            {
                axis_center
            } else {
                (axis_range.0 + axis_range.1) * 0.5
            };

            let center_ratio = (center_used - axis_range.0) / range;
            let new_min = center_used - new_range * center_ratio;
            let new_max = center_used + new_range * (1.0 - center_ratio);

            if new_max > new_min {
                Some((new_min, new_max))
            } else {
                None
            }
        };

        if let Some(new_x) = zoom_axis(self.view_range, center[0]) {
            self.view_range = new_x;
        }
        if let Some(new_y) = zoom_axis(self.view_y_range, center[1]) {
            self.view_y_range = new_y;
        }
    }

    /// Масштабирует отображение под выделенную область
    pub fn fit_to_selection(&mut self) {
        // Проверяем, есть ли активное выделение
        if let Some(rect) = &self.selection_rect {
            // Нормализуем прямоугольник выделения
            let normalized = rect.normalized();
            // Устанавливаем диапазон по X равным границам выделения с небольшим отступом
            let x_min = normalized.start[0];
            let x_max = normalized.end[0];
            if (x_max - x_min).abs() > 1e-9 {
                // Вычисляем отступ (5% от ширины выделения)
                let margin = (x_max - x_min) * 0.05;
                let new_min = x_min - margin;
                let new_max = x_max + margin;
                // Применяем новые границы, если они корректны
                if new_max > new_min {
                    self.view_range.0 = new_min;
                    self.view_range.1 = new_max;
                }
            }
            // Также устанавливаем диапазон по Y равным границам выделения с отступом
            let y_min = normalized.start[1];
            let y_max = normalized.end[1];
            if (y_max - y_min).abs() > 1e-9 {
                // Вычисляем отступ по Y (5% от высоты выделения)
                let margin_y = (y_max - y_min) * 0.05;
                let new_y_min = y_min - margin_y;
                let new_y_max = y_max + margin_y;
                // Применяем новые границы по Y, если они корректны
                if new_y_max > new_y_min {
                    self.view_y_range.0 = new_y_min;
                    self.view_y_range.1 = new_y_max;
                }
            }
        }
    }
}
