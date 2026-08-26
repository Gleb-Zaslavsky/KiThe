# Codex Task: заменить phase-stability subsystem KiThe на общую TPD-постановку

## Цель

Полностью заменить текущий ограниченный механизм phase stability/activation на математически корректный tangent-plane stability analysis для многокомпонентных фаз.

Реализация должна соответствовать документу **«KiThe: математический и архитектурный контракт фазовой устойчивости»**.

Это рефакторинг с изменением архитектуры, а не задача обратной совместимости.

---

## 1. Правило работы со старым кодом

Существующий phase-stability код НЕ является API-контрактом.

Разрешается и требуется:

* удалять устаревшие функции;
* удалять устаревшие enum variants;
* менять структуры и их поля;
* менять сигнатуры внутренних функций;
* переписывать call sites;
* удалять tests, проверяющие старую неправильную семантику, и заменять их тестами новой семантики.

НЕ делать:

* compatibility adapters;
* legacy wrappers;
* deprecated copies старых функций;
* параллельные `old/new` реализации;
* feature flags для сохранения старого stability criterion;
* traits или abstraction layers только ради возможности вызвать старый код;
* fallback к старому pure-phase criterion при ошибке нового алгоритма.

После рефакторинга должна существовать одна актуальная реализация phase stability.

Если старый код противоречит новому mathematical contract — удалить его.

---

## 2. Сначала провести dependency audit

До изменения кода найти все определения и call sites, связанные как минимум с:

```text
compute_phase_stability_reports
PhaseStabilityReport
PhaseStabilityModel
ElementPotentialReport
driving_force
dg_create
dg_keep
seed_activated_phase
PhaseSeedPolicy
PhaseManager
PreparedActiveSetCandidate
PreparedPhaseControlRunner
P,H active-set solve
complementarity reports
phase-control tests
```

Также найти код сериализации/экспорта/reporting, если он зависит от старых полей.

Составить внутреннюю карту зависимостей и затем менять цепочку целиком.

Не сохранять старую сигнатуру функции только потому, что у нее много call sites.

---

## 3. Выделить stability analysis в отдельный модуль

Создать отдельный модуль, например

```text
equilibrium_phase_stability.rs
```

и перенести туда физику фазовой устойчивости.

`equilibrium_workflows.rs` должен отвечать за внешний phase-control workflow, а не содержать реализацию TPD mathematics.

В новом модуле должны находиться:

* reconstruction of elemental potentials;
* construction of candidate-phase TPD problem;
* elemental-direction feasibility;
* minimization of TPD;
* stability report structures.

Не переносить туда:

* изменение `PhaseSet`;
* hysteresis state machine;
* cycle detection;
* transactional publish.

---

## 4. Удалить старую физическую классификацию

Текущий подход вида

```text
FixedIdealGas
PureCondensedSpecies
Unsupported
```

не должен задавать фундаментальную структуру stability analysis.

Удалить `PhaseStabilityModel`, если после рефакторинга он больше не несет самостоятельного смысла.

Использовать существующую thermodynamic phase/activity model.

Чистая однокомпонентная конденсированная фаза должна автоматически работать как случай `q = 1` общей TPD-задачи.

Не оставлять отдельный pure-phase algorithm рядом с multicomponent algorithm.

---

## 5. Вычислять химические потенциалы из канонической модели

Для принятого состояния вычислять

```text
μ_i = g_i°(T) + RT ln a_i
```

через тот же activity-model boundary, который используется основной equilibrium formulation.

Не дублировать формулы `IdealGas`/`IdealSolution` внутри stability subsystem, если каноническая функция activity evaluation уже существует.

Если для корректной TPD-минимизации нужен новый метод activity model, добавить его в канонический `PhaseActivityModel`, а не создавать вторую независимую реализацию термодинамики.

---

## 6. Реконструировать элементные потенциалы

Для активного принятого состояния решить

```text
A_active λ ≈ μ_active
```

устойчивым методом на основе SVD или существующей надежной LA-инфраструктуры.

Сохранить diagnostics:

```text
potentials
rank
max_abs_residual
```

Проверять, что residual соответствует acceptance tolerance.

ВАЖНО:

`λ` может быть неединственным при rank-deficient `A_active`.

Не считать minimum-norm `λ` физически достаточным для произвольного candidate elemental vector.

---

## 7. Ввести elemental-direction feasibility

Для каждой candidate phase и trial composition `x` вычислять

```text
c(x) = A_phase^T x
```

Пробное появление допустимо только если

```text
c(x) ∈ Range(A_active^T)
```

с заданным численным tolerance.

Реализовать это через rank/null-space/SVD formulation.

Предпочтительная форма:

```text
N^T c(x) = 0
```

где столбцы `N` задают `Null(A_active)`.

Если active elemental rank полный, это ограничение автоматически исчезает.

Не отбрасывать rank-deficient Active Set как ошибочный только ради упрощения stability calculation.

Не вычислять TPD для physically infeasible trial composition.

---

## 8. Реализовать общий TPD objective

Для candidate phase `β` и состава `x`:

```text
TPDβ(x)
    = Σ_i x_i [ μ_i^β(x,T,P) - a_i·λ ]
```

при

```text
x_i >= 0
Σ_i x_i = 1
elemental-direction feasibility
```

Все величины TPD хранить в J/mol.

Не вводить arbitrary trace amount в definition of TPD.

TPD — молярный differential criterion, а не finite-seed experiment.

---

## 9. Для текущих IdealSolution/IdealGas использовать точную структуру задачи

Для `IdealSolution`:

```text
μ_i = g_i° + RT ln(x_i)
```

Для `IdealGas`:

```text
μ_i = g_i° + RT ln(x_i P/P0)
```

Если elemental rank полный, допускается и предпочтительна аналитическая/специализированная реализация глобального минимума вместо generic optimizer.

Если присутствуют дополнительные linear feasibility constraints из rank-deficient active assemblage, решить соответствующую constrained convex problem корректно.

Не использовать:

* penalty-only constraints;
* normalization after an unconstrained optimum, если это нарушает elemental feasibility;
* случайный composition seed как замену minimization.

Если в RustedSciThe уже существует подходящий constrained optimizer/root solver — использовать его.

Если его нет, не ослаблять математическую постановку. Реализовать минимально необходимый численный метод в том слое, где ему архитектурно место. Не добавлять стороннюю numerical dependency без необходимости.

---

## 10. Подготовить архитектуру к non-ideal models

Новый stability subsystem должен формулировать задачу через:

```text
phase composition
      ↓
activity model
      ↓
chemical potentials
      ↓
TPD
```

Не зашивать assumption `ln(a_i)=ln(x_i)` в общий слой.

Для будущей non-ideal model TPD может стать non-convex и требовать multistart/global strategy.

Не реализовывать сейчас фиктивную non-ideal модель.

Но структура данных и objective boundary не должны требовать нового архитектурного рефакторинга при ее появлении.

---

## 11. Новый PhaseStabilityReport

Перепроектировать отчет примерно вокруг следующей семантики:

```rust
PhaseStabilityReport {
    phase,
    active,
    minimum_tpd,
    incipient_composition,
    element_potentials,
    ...
}
```

Точные Rust-типы выбрать по существующим typed IDs и style проекта.

Обязательные свойства:

`minimum_tpd`

* `Some(value)` для реально выполненного stability test;
* J/mol.

`incipient_composition`

* соответствует порядку species внутри candidate phase;
* нормирована;
* является `argmin TPD`;
* для pure phase равна `[1.0]`.

Не хранить поле с названием `driving_force`, если фактически теперь вычисляется global minimum TPD. Использовать терминологию, отражающую математику.

Не сохранять `PhaseStabilityModel::PureCondensedSpecies` ради старых call sites.

---

## 12. Использовать minimizer composition для phase seed

Изменить phase activation pipeline.

Сейчас phase seed не должен распределяться между компонентами новой фазы произвольно.

Новый dataflow:

```text
stability analysis
       │
       ├── minimum_tpd
       │       ↓
       │   activation decision
       │
       └── x*
               ↓
        phase seed composition
```

`PhaseSeedPolicy` продолжает определять общий seed amount:

```text
n_phase_seed
```

Но species amounts задаются:

```text
n_i_seed = n_phase_seed * x_i*
```

Для координат `x_i*=0`, несовместимых со строгой log-mole positivity, применить очень малый composition floor и перенормировать.

Не заменять `x*` равномерным распределением.

Изменить `seed_activated_phase` или заменить его новой функцией с корректной семантикой.

Если старая функция после этого не нужна — удалить.

---

## 13. Интеграция с hysteresis

Существующие thresholds:

```text
dg_create
dg_keep
```

не являются определением stability criterion.

Они должны применяться к `minimum_tpd` после его вычисления.

Сохранить общий смысл:

```text
minimum_tpd < dg_create
```

→ inactive phase is sufficiently unstable to be activated.

Не менять hysteresis architecture больше, чем необходимо для перехода от старого `driving_force` к `minimum_tpd`.

Более глубокий hysteresis refactor будет отдельной задачей.

---

## 14. P,H integration

Для monolithic `P,H=const` stability analysis выполнять при:

```text
T = accepted_solution.temperature
P = problem.pressure
```

Не использовать initial temperature seed.

Если candidate phase активируется:

1. использовать `incipient_composition` для phase seed;
2. сохранить accepted temperature как continuation temperature seed;
3. перестроить Active Set;
4. запустить новый monolithic P,H solve;
5. позволить ему заново определить и composition, и temperature при прежнем `H_target`.

Не создавать отдельный P,H stability criterion.

---

## 15. Excluded phases

`Excluded` phase не должна участвовать в TPD minimization.

Не создавать для нее fake report с `Unsupported`.

Либо не включать ее в список stability calculations, либо выдавать отчет со строго определенным `NotEvaluated/Excluded` состоянием, только если reporting architecture действительно этого требует.

Не путать:

```text
excluded
```

и

```text
thermodynamic model unsupported.
```

Это разные причины отсутствия результата.

---

## 16. Tests: обязательный минимум

Добавить/переписать unit и integration tests.

### A. Pure phase regression

Однокомпонентная ideal condensed phase.

Проверить, что общий TPD algorithm дает тот же результат, что аналитическая формула:

```text
g° - a·λ
```

Отдельный pure-phase implementation при этом отсутствует.

### B. Stable multicomponent ideal solution

Двух- или трехкомпонентная inactive phase.

Проверить:

```text
minimum_tpd > 0
```

и отсутствие activation.

### C. Unstable multicomponent ideal solution

Сконструировать систему с известным отрицательным минимумом.

Проверить:

```text
minimum_tpd < 0
```

и корректную `incipient_composition`.

### D. Known analytical ideal-solution minimum

Для full-rank случая сравнить numerical result с аналитическим минимумом ideal mixture.

Проверить и `TPD_min`, и `x*`.

### E. Phase seed from minimizer

Проверить, что после activation species seed ratios соответствуют `x*`, а не равномерному распределению.

### F. Boundary composition

Случай, когда optimum находится близко к границе simplex.

Проверить:

* корректный `TPD_min`;
* отсутствие `ln(0)`;
* корректный composition floor только на стадии seed.

### G. Rank-deficient active elemental system

Построить Active Set с dependent element columns.

Проверить, что:

* λ reconstruction допускает nonunique potentials;
* feasible candidate direction дает invariant TPD;
* результат не зависит от arbitrary null-space component λ.

### H. Elementally infeasible candidate direction

Проверить, что trial composition вне `Range(A_active^T)` не получает фиктивную driving force.

### I. P,H

Проверить, что stability calculation использует solved temperature, а не initial seed.

### J. Hysteresis integration

Проверить существующие create/keep threshold transitions уже на `minimum_tpd`.

### K. Existing phase-control regressions

Все существующие physical regression tests P,T и P,H должны пройти.

Тесты, фиксирующие старое отсутствие поддержки multicomponent phase stability, удалить как устаревшие.

---

## 17. Проверки инвариантов

После каждого stability solve проверять:

```text
Σ x_i = 1
x_i >= 0
minimum_tpd finite
incipient composition finite
elemental feasibility residual within tolerance
element-potential reconstruction residual within tolerance
```

Для phase seed дополнительно:

```text
all seeded n_i > 0
Σ n_i = requested phase seed amount
```

с учетом численного tolerance.

---

## 18. Не делать скрытых fallback

Если новый TPD solver не сошелся:

НЕ:

* использовать старый pure-phase test;
* брать равномольный trial composition;
* считать phase stable;
* считать phase unstable;
* активировать фазу «на всякий случай».

Вернуть явную диагностируемую ошибку stability calculation.

Failure to prove stability is not proof of stability.

---

## 19. Удаление obsolete code

После успешной миграции:

* удалить старую branch `phase.species.len() != 1 => ValidationNotApplicable`;
* удалить специальный pure condensed stability path;
* удалить obsolete `PhaseStabilityModel` variants/enum, если он более не нужен;
* удалить старые helper functions;
* удалить compatibility wrappers;
* удалить dead imports;
* удалить tests старой семантики;
* выполнить `cargo fmt`;
* выполнить `cargo clippy`;
* выполнить полный relevant test suite.

Не оставлять комментарии вида:

```rust
// old implementation kept for compatibility
```

если compatibility явно не требуется внешним public API.

Если public API действительно затронут, изменить его осознанно и обновить call sites/documentation. Не строить внутренний adapter только для сохранения ошибочной модели.

---

## 20. Критерий завершения

Задача выполнена только если dataflow имеет вид:

```text
accepted fixed-set solution
           │
           ▼
chemical potentials
           │
           ▼
element-potential reconstruction
           │
           ▼
candidate phase
           │
           ▼
TPD minimization over
admissible composition
           │
       ┌───┴────┐
       ▼        ▼
   TPD_min      x*
       │        │
       ▼        ▼
 hysteresis   phase seed
       │        │
       └───┬────┘
           ▼
      new Active Set
           │
           ▼
     new fixed-set solve
```

Многокомпонентная `IdealSolution` должна участвовать в phase creation на тех же основаниях, что и однокомпонентная фаза.

Никакой временной ветки «multicomponent stability unsupported» после завершения задачи оставаться не должно.
