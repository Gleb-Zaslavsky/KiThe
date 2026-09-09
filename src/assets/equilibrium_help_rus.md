# Калькулятор химического равновесия

Форма создает один канонический запрос `EquilibriumCalculator`; отдельного
GUI-солвера нет. Результаты являются неизменяемыми принятыми snapshots.

## Setup

`P,T = const` ищет состав при фиксированных давлении и температуре. `P,H =
const` ищет состав и температуру при фиксированных давлении и полной
энтальпии. `Point` решает одну точку, `Range` доступен для P,T, сохраняет
направление сетки и переиспользует formulation и последнее принятое состояние.

`Pressure [Pa]` -- давление системы, конечное и положительное. `Reference
pressure [Pa]` -- standard-state pressure в законе активности идеального газа; оно
не заменяется автоматически давлением системы и должно соответствовать
библиотеке. `Temperature [K]` -- температура задачи P,T.

Для P,T range поля `Start`, `End`, `Points` задают сетку; нужно минимум две
точки. В P,H `Target total enthalpy [J]` -- extensive enthalpy всего inventory,
не J/mol. `Lower bound`, `Upper bound`, `Initial seed [K]` задают интервал и
начальную оценку поиска температуры; seed должен быть внутри интервала, а
интервал -- внутри области термохимических данных.

`Explicit species` принимает точные записи. `Search by elements` ищет кандидатов,
но их еще нужно назначить фазе. `Elements` задает набор элементов; `Exact`
требует тот же набор, `Subset` разрешает подмножество. Фильтры состояния и
лимит кандидатов делают поиск в каталоге детерминированным.

Оба режима проходят через публичный типизированный facade
`ChemEquilibrium::prelude`. Explicit species задает молекулярное начальное
состояние; element mode задает сохраняемый элементный inventory, а поиск
кандидатов только выбирает реальные фазовые вещества.

Для каждой фазы `Phase` является стабильным identity, `Physical state` задает
gas/liquid/solid/condensed semantics, а `Model` выбирает ideal gas,
ideal solution или pure condensed. Несовместимые пары отклоняются. `Substance`
-- идентификатор библиотеки, `Initial moles` -- количество именно этой фазы.
`Add phase`, `Add component`, `Remove` меняют документ и требуют новой
подготовки.

`Validate document` выполняет field-local проверку. `Prepare canonical request`
передает validated state фасаду, но не запускает solve. `Run prepared request`
запускает worker. `Cancel` не позволяет позднему или устаревшему выводу заменить
принятый результат.

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Задача, точка/диапазон, давление и температура | Положительные Pa и K; по умолчанию `Point` | Всегда; диапазон только для P,T | Инвалидирует запрос |
| Целевая энтальпия P,H и температурный интервал | Экстенсивные J; конечный интервал K и seed внутри него | Только P,H | Инвалидирует запрос |
| Вещества, элементы, фазы и начальные моли | Точные идентификаторы; конечные неотрицательные mol | Согласно режиму задания | Инвалидирует запрос |
| Validate / Prepare / Run / Cancel | Сначала проверка; подготовка создаёт неизменяемый запрос | Доступность зависит от документа и worker | Подготовка / выполнение |

## Phase control

`Fixed declared phases` сохраняет ровно объявленный набор фаз. `Bounded phase
control` включает outer lifecycle, который может создавать или удалять
кандидаты после TPD/stability checks.

`Phase epsilon [mol]` -- threshold присутствия фазы. `Creation driving force`
управляет появлением отсутствующей фазы, `Keep driving force` -- сохранением
активной. Их разность образует hysteresis и предотвращает chatter. `Maximum
phase iterations` ограничивает работу; непроверенное промежуточное состояние
никогда не публикуется.

`Positive inventory` начинает с фаз с положительным количеством. `All declared
candidates` проверяет также declared phases с нулевым количеством, что удобно
для тестов появления. Это lifecycle policy, а не новая физическая модель.

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Fixed / bounded lifecycle | По умолчанию фиксированный набор объявленных фаз | Всегда | Инвалидирует запрос |
| Phase epsilon | Положительный порог mol | Bounded lifecycle | Инвалидирует запрос |
| Create / keep driving force | Конечные пороги J/mol; разность образует hysteresis | Bounded lifecycle | Инвалидирует запрос |
| Initial phase set и лимит итераций | Positive inventory; конечное положительное число итераций | Bounded lifecycle | Инвалидирует запрос |

## Libraries

`Engine default` передает priority общему repository. `Explicit policy`
использует упорядоченный `Priority libraries` и закрытый `Permitted libraries`.
`Offline mode` запрещает сеть. `Allow online NIST fallback` является явным
opt-in и не подменяет локальную запись молча.

`Library catalog status` показывает loading, loaded, unavailable или failed для
общего каталога. `Lookup instruction` описывает желаемый источник до resolve.
Фактические record identity, library, version и provenance находятся в
принятых Results.

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Lookup policy | По умолчанию policy движка, пока не выбран explicit режим | Всегда | Инвалидирует запрос |
| Priority и permitted libraries | Упорядоченный приоритет и закрытый allow-list | Explicit policy | Инвалидирует запрос |
| Offline / NIST fallback | Offline безопасен по умолчанию; NIST требует opt-in | Всегда | Инвалидирует запрос |
| Load catalog | Действие над общим локальным repository | Пока каталог не готов | Только подготовка |

## Numerics

`Production default cascade` использует production backend cascade. `Single backend`
удобен для диагностики и не имеет fallback. `Custom cascade` предназначен для
контролируемых сравнений и должен содержать уникальные поддерживаемые backends.

`P,H route` предлагает `Auto`, `Monolithic` и `Nested temperature`, если маршрут
поддержан. `Tolerance` определяет численное принятие, `Maximum iterations`
ограничивает работу. Пустые overrides используют production policy. `Enable
scaling` меняет численные координаты, но не физические единицы. `Trace seed
policy` задает положительные log-mole seeds и не добавляет inventory.
`Cascade total iterations` ограничивает общую работу fallback, сохраняя evidence.

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Solver policy | По умолчанию production cascade | Всегда | Инвалидирует запрос |
| Backend и custom cascade | Поддерживаемые уникальные backends | Single backend / custom cascade | Инвалидирует запрос |
| P,H route | По умолчанию Auto | Только P,H | Инвалидирует запрос |
| Tolerance, iteration и scaling overrides | Пустое поле использует production defaults; overrides конечны и положительны | Advanced numerics | Инвалидирует запрос |
| Trace-seed policy | Не добавляет физический inventory | Advanced numerics | Инвалидирует запрос |

## Output

`Collect timing` записывает timing этапов и точек. `Phase lifecycle trace` хранит
TPD, stability, activation, deactivation и transitions bounded control.
`Maximum retained lifecycle events` ограничивает память trace, сохраняя последние
события. `Range lifecycle trace` выбирает endpoints, transitions или every point.

`K_eq validation` включает независимую проверку малых систем и не заменяет
общий production solver. `Result basis` выбирает moles, mole fractions или
phase totals для display series.

`Result table` выбирает detailed moles + fractions или compact moles. `Hide below
mole fraction` фильтрует строки видимой таблицы; ноль показывает все компоненты.
Обе настройки только display-only и не меняют snapshot.

`Plot target` выбирает embedded plot, KiThePlot, оба варианта или none. `Y scale`
выбирает linear или log10. `PCHIP display resampling` доступен для P,T range;
число точек, interpolation space и clamp меняют только display data, а не точки
решения.

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Timing и lifecycle trace | По умолчанию выключены; лимит событий ограничивает память | Diagnostics options | Диагностика выполнения |
| K_eq validation | По умолчанию выключена; применима только к малым системам | Diagnostics options | Инвалидирует запрос |
| Таблица, basis, порог и видимые series | Параметры представления; конечный порог доли | Accepted results | Только отображение |
| Plot target и scale | Без графика / linear по умолчанию | Accepted results | Только отображение |
| PCHIP resampling | Выключен; число display points для диапазона P,T | Results P,T range | Только отображение |

## Results

### Контракт элементов управления

| Элемент | Область / эффективное значение | Когда виден | Класс изменения |
| --- | --- | --- | --- |
| Accepted snapshot | Неизменяемое принятое физическое состояние | После успешного solve | Только отображение |
| Diagnostics, provenance и lifecycle evidence | Неизменяемые доказательства, сохранённые запросом | Если включён соответствующий сбор | Только отображение |
| Range selection и plots | Производны только от принятых точек диапазона | Принятый P,T range | Только отображение |

Каждая принятая точка является immutable accepted snapshot и доступна только для чтения. `Solve diagnostics` содержит residuals,
element conservation, backend attempts, acceptance, fallback reasons, timing и
validation evidence. `Lookup provenance` показывает resolved source каждого
компонента. `Phase lifecycle` показывает TPD и transitions, если trace сохранен.
Range summary содержит direction, accepted points, formulation builds/reuses,
normalization recoveries, transitions и timing. Неудачные точки не выдаются за
принятые, а plotting не запускает solver повторно.
