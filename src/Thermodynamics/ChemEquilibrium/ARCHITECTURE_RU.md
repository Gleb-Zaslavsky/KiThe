# Архитектурный обзор модуля ChemEquilibrium

## Содержание

1. [Введение и общая архитектура](#1-введение-и-общая-архитектура)
2. [Основные структуры данных](#2-основные-структуры-данных)
3. [Dataflow: путь данных от входа до решения](#3-dataflow-путь-данных-от-входа-до-решения)
4. [Архитектурно-математические особенности](#4-архитектурно-математические-особенности)
   - [4.1 Log-mole формулировка](#41-log-mole-формулировка)
   - [4.2 Символьный и численный стек](#42-символьный-и-численный-стек)
   - [4.3 Масштабирование (Row scaling и Variable scaling)](#43-масштабирование-row-scaling-и-variable-scaling)
   - [4.4 Предотвращение нефизических значений](#44-предотвращение-нефизических-значений)
   - [4.5 SVD реакционный базис](#45-svd-реакционный-базис)
   - [4.6 Управление фазами: гистерезис, No flip-flop, bounded outer loop](#46-управление-фазами-гистерезис-no-flip-flop-bounded-outer-loop)
   - [4.7 Phase stability (проверка стабильности фаз)](#47-phase-stability-проверка-стабильности-фаз)
   - [4.8 Element-validation after activation](#48-element-validation-after-activation)
   - [4.9 Независимая K_eq валидация](#49-независимая-k_eq-валидация)
   - [4.10 Транзакционность публикации состояния](#410-транзакционность-публикации-состояния)
   - [4.11 Политика решателя (SolverPolicy) и каскад](#411-политика-решателя-solverpolicy-и-каскад)
   - [4.12 Температурные серии и постобработка](#412-температурные-серии-и-постобработка)
5. [Навигационная карта](#5-навигационная-карта)
   - [5.1 С чего начать чтение](#51-с-чего-начать-чтение)
   - [5.2 Карта файлов](#52-карта-файлов)
   - [5.3 Типичные сценарии использования](#53-типичные-сценарии-использования)

---

## 1. Введение и общая архитектура

Модуль `ChemEquilibrium` — это ядро расчёта химического равновесия в библиотеке KiThe. Он решает задачу минимизации энергии Гиббса для многофазных систем при фиксированных температуре и давлении.

Архитектура модуля построена вокруг **трёх ключевых слоёв**:

```
┌─────────────────────────────────────────────────────────────┐
│                    Публичный фасад                           │
│  phase_equilibrium_workflow.rs, easy_equilibrium.rs          │
│  gas_solver(), solve_resolved_pt()                           │
├─────────────────────────────────────────────────────────────┤
│                 Канонический слой                            │
│  equilibrium_problem.rs → EquilibriumProblem                │
│  equilibrium_log_moles.rs → EquilibriumLogMoles             │
│  equilibrium_workflows.rs → PhaseManager                    │
├─────────────────────────────────────────────────────────────┤
│              Слой бэкенда (адаптеры)                         │
│  equilibrium_rst_backend.rs   (RustedSciThe — символьный)   │
│  equilibrium_legacy_backend.rs (Legacy — численный)          │
│  equilibrium_backend_adapter.rs (единый интерфейс)           │
├─────────────────────────────────────────────────────────────┤
│              Слой валидации                                  │
│  equilibrium_validation.rs   (каноническая проверка)         │
│  equilibrium_constant_solver.rs (независимый K_eq решатель)  │
│  equilibrium_constant_cross_validation.rs (сравнение)        │
├─────────────────────────────────────────────────────────────┤
│              Мостовой слой (фазовая подсистема)              │
│  phase_equilibrium_problem.rs → мост ResolvedPhaseSystem    │
│  phase_equilibrium_solution.rs → MultiphaseEquilibriumSolution│
│  equilibrium_multiphase_domain.rs → типобезопасная граница   │
└─────────────────────────────────────────────────────────────┘
```

**Ключевые принципы архитектуры:**

1. **Типобезопасность** — идентификаторы SpeciesId, ElementId, PhaseIndex, ReactionId исключают путаницу между индексами разных сущностей.
2. **Явная политика решателя** — выбор бэкенда (RST или Legacy) и порядок каскада задаются типизированной политикой, а не цепочкой if-else.
3. **Транзакционность** — ни один бэкенд не публикует состояние напрямую; решение проходит валидацию и только затем публикуется.
4. **Независимая валидация** — K_eq решатель использует отдельную формулировку (закон действующих масс), не пересекающуюся с основной residual/Jacobian.
5. **Управление фазами** — bounded outer loop с гистерезисом, детекцией циклов и typed transition plan.

---

## 2. Основные структуры данных

### 2.1 Типобезопасные идентификаторы

Файл: [`equilibrium_ids.rs`](equilibrium_ids.rs)

Четыре обёртки над `usize`, каждая со своей семантикой:

| Тип | Назначение |
|-----|------------|
| `SpeciesId` | Индекс вещества в упорядоченном списке решателя |
| `ElementId` | Индекс химического элемента |
| `PhaseIndex` | Плотный индекс фазы в массивах равновесия |
| `ReactionId` | Индекс независимой реакции |

Все создаются через `::new(index, upper_bound)` с проверкой границ.

### 2.2 Граница задачи (Problem boundary)

Файл: [`equilibrium_problem.rs`](equilibrium_problem.rs)

```rust
EquilibriumConditions {
    temperature: f64,       // K
    pressure: f64,          // Pa
    reference_pressure: f64,// Pa (обычно 101325)
}

EquilibriumProblem {
    components: Vec<EquilibriumComponentDescriptor>, // фазово-квалифицированные компоненты
    labels: Vec<String>,                             // метки для отчётов
    initial_moles: Vec<f64>,                         // физические начальные моли
    initial_log_moles: LogMolesInitialGuess,         // ln(n_i) начальное приближение
    element_composition: DMatrix<f64>,               // матрица элементного состава (вещества × элементы)
    gibbs: Vec<GibbsFn>,                             // функции G0(T) для каждого вещества
    phases: Vec<Phase>,                              // описание фаз
    conditions: EquilibriumConditions,               // T, P, P0
}

PreparedEquilibriumProblem {
    problem: EquilibriumProblem,
    reaction_basis: ReactionBasis,     // SVD-базис независимых реакций
    element_totals: Vec<f64>,          // полные моли каждого элемента
    species_phase: Vec<usize>,         // отображение вещество → фаза
    phase_stoichiometry: Vec<Vec<f64>>,// стехиометрия фаз
}
```

`EquilibriumProblem` — это **входная граница**: все поля проверены на конечность и физичность. `PreparedEquilibriumProblem` — это подготовленная форма, содержащая все матрицы и порядки, необходимые любому бэкенду.

### 2.3 Решение (Solution)

```rust
EquilibriumSolution {
    log_moles: Vec<f64>,           // ln(n_i) — невязка решателя
    moles: Vec<f64>,               // физические моли
    conditions: EquilibriumConditions,
    validation: EquilibriumCandidateReport,  // отчёт о проверке
}
```

`EquilibriumSolution` — **иммутабельный результат**. Он содержит не только численное решение, но и доказательство его приёмки (validation report).

### 2.4 Основной оркестратор

Файл: [`equilibrium_log_moles.rs`](equilibrium_log_moles.rs) (3881 строка)

```rust
EquilibriumLogMoles {
    subs_data: SubsData,                    // термохимические данные
    elem_composition: DMatrix<f64>,         // элементный состав
    reaction_basis: ReactionBasis,          // SVD-базис
    initial_guess: Option<Vec<f64>>,        // начальное приближение (ln)
    P: f64, T: f64, p0: f64,               // условия
    stoich_matrix: DMatrix<f64>,            // матрица стехиометрии
    n0: Vec<f64>,                           // начальные моли
    gibbs: Vec<GibbsFn>,                    // G0(T) замыкания
    gibbs_sym: Vec<Expr>,                   // символьные G0(T)
    phases: Vec<Phase>,                     // фазы
    species_phase: Vec<usize>,              // вещество → фаза
    phase_active_mask: Vec<bool>,           // маска активных фаз
    solver_settings: EquilibriumSolverSettings,
    phase_manager: PhaseManager,
    // ... поля для температурных серий, отчётов, кэшей
}
```

Это **исторический фасад**, который постепенно вытесняется типобезопасными API. Содержит как публичные методы (`solve`, `solve_for_T_range`), так и внутренние (`solve_candidate_from_seed`, `solve_fixed_active_set_candidate`).

### 2.5 Настройки решателя

```rust
EquilibriumSolverSettings {
    solver_params: SolverParams,         // max_iter, tol, lambda, alpha_min, delta...
    solver: Solvers,                     // LM, NR, TR (legacy)
    scaling_flag: bool,                  // включить row scaling
    trace_seed_policy: TraceSpeciesSeedPolicy,  // политика начального приближения
    solver_policy: Option<SolverPolicy>,        // Single/Cascade
    solver_budget: Option<SolverCascadeBudget>, // бюджет каскада
    keq_validation_mode: EquilibriumConstantValidationMode,
    keq_validation_tolerances: EquilibriumConstantCrossValidationTolerances,
    continuation_seed_policy: ContinuationSeedPolicy,
}
```

### 2.6 Политика решателя

Файл: [`equilibrium_solver_policy.rs`](equilibrium_solver_policy.rs)

```rust
enum SolverBackend {
    Legacy(Solvers),           // LM, NR, TR
    RustedSciThe(RustedSciTheSolver),  // 6 методов
}

enum SolverPolicy {
    Single(SolverBackend),     // ровно один бэкенд
    Cascade(Vec<SolverBackend>), // упорядоченный каскад
}

SolverCascadeBudget {
    max_attempts: usize,
    max_iterations_per_attempt: usize,
    max_total_iterations: usize,
}
```

### 2.7 Управление фазами

Файл: [`equilibrium_workflows.rs`](equilibrium_workflows.rs)

```rust
PhaseManager {
    phase_eps: f64,                          // порог уничтожения фазы
    phase_hysteresis: PhaseHysteresisPolicy, // гистерезис
    max_phase_iterations: usize,             // макс. итераций outer loop
    initial_phase_set: InitialPhaseSet,      // политика начального набора
}

PhaseSet {
    statuses: Vec<PhaseStatus>,  // Active, Inactive, Excluded, Appeared, Disappeared
}

PhaseTransitionPlan {
    Deactivate { phase: PhaseIndex },
    Activate { phase: PhaseIndex },
    Hold { phase: PhaseIndex },
    NoTransition { reason: PhaseStabilityReason },
}
```

### 2.8 Проекция активного набора

Файл: [`equilibrium_active_set.rs`](equilibrium_active_set.rs)

```rust
ActiveSetProjection {
    active_phases: Vec<PhaseIndex>,
    active_species: Vec<SpeciesId>,
    phases: Vec<Phase>,
    species_phase: Vec<usize>,
    element_composition: DMatrix<f64>,
    reaction_basis: ReactionBasis,
    full_species_count: usize,
}
```

Этот тип **замораживает** решение outer loop для одного каскада бэкенда. Он владеет единственной картой отображения глобальных индексов в локальные.

### 2.9 Мостовые структуры

Файл: [`phase_equilibrium_problem.rs`](phase_equilibrium_problem.rs)

```rust
PhaseEquilibriumMetadata {
    layout: SystemLayout,
    components: Vec<EquilibriumComponentDescriptor>,
    phases: Vec<EquilibriumPhaseDescriptor>,
    provenance: ResolvedPhaseSystemReport,
}

PhaseEquilibriumBuildRequest {
    resolved: &ResolvedPhaseSystem,
    conditions: EquilibriumConditions,
    initial_composition: MultiphaseInitialComposition,
    trace_seed_policy: TraceSpeciesSeedPolicy,
    model_policy: SupportedPhaseModelPolicy,
}
```

Файл: [`phase_equilibrium_solution.rs`](phase_equilibrium_solution.rs)

```rust
MultiphaseEquilibriumSolution {
    metadata: PhaseEquilibriumMetadata,
    build_report: PhaseEquilibriumBuildReport,
    accepted_solution: EquilibriumSolution,
    phase_totals: Vec<f64>,
    mole_fractions: Vec<f64>,
    phase_statuses: Vec<PhaseStatus>,
    solve_report: EquilibriumSolveReport,
    keq_validation_status: Option<EquilibriumConstantCrossValidationStatus>,
}
```

### 2.10 Многофазная доменная граница

Файл: [`equilibrium_multiphase_domain.rs`](equilibrium_multiphase_domain.rs)

```rust
MultiphaseEquilibriumLayout {
    phase_specs: Vec<PhaseSpec>,
    system_layout: SystemLayout,
    fingerprint: u64,
}

MultiphaseInitialComposition {
    layout_fingerprint: u64,
    moles: Vec<f64>,
}
```

### 2.11 Валидация

Файл: [`equilibrium_validation.rs`](equilibrium_validation.rs)

```rust
EquilibriumCandidateReport {
    residual_l2_norm: f64,
    residual_rms: f64,
    max_abs_residual: f64,
    raw_residual_l2_norm: f64,
    max_abs_element_balance_error: f64,
    reaction_affinity_l2_norm: f64,
    min_moles: f64,
    // ...
}

EquilibriumAcceptanceCriteria {
    residual_tolerance: f64,
    element_balance_tolerance: f64,
    reaction_affinity_tolerance: f64,
}
```

### 2.12 Независимый K_eq решатель

Файл: [`equilibrium_constant_solver.rs`](equilibrium_constant_solver.rs)

```rust
EquilibriumConstantSolver {
    settings: EquilibriumConstantSolverSettings,
    validation_tolerances: EquilibriumConstantValidationTolerances,
    mode: EquilibriumConstantSolverMode,  // Off, WhenApplicable, Required
}

EquilibriumConstantSolveResult {
    extent: f64,
    moles: Vec<f64>,
    validation: EquilibriumConstantValidationReport,
    report: EquilibriumConstantSolveReport,
}
```

---

## 3. Dataflow: путь данных от входа до решения

### 3.1 Полный путь для многофазной системы (рекомендуемый)

```
ResolvedPhaseSystem (фазовая подсистема)
    │  [phase_equilibrium_problem.rs]
    ▼
PhaseEquilibriumBuildRequest
    │  resolved, conditions, initial_composition,
    │  trace_seed_policy, model_policy
    │  [phase_equilibrium_problem.rs → build_phase_equilibrium_problem()]
    ▼
PhaseEquilibriumMetadata  ←─── MultiphaseEquilibriumLayout
    │  (валидация: дубликаты фаз, поддержка моделей,
    │   сортировка по PhaseId, fingerprint)
    │  [equilibrium_multiphase_domain.rs]
    ▼
EquilibriumProblem
    │  (проверка: конечность, размерности, seed policy)
    │  [equilibrium_problem.rs → EquilibriumProblem::validate()]
    ▼
PreparedEquilibriumProblem
    │  (SVD реакционного базиса, element totals,
    │   species_phase map, phase_stoichiometry)
    │  [equilibrium_problem.rs → PreparedEquilibriumProblem::new()]
    ▼
EquilibriumSolveContract
    │  (замыкания: residual, jacobian, rst_problem)
    │  [equilibrium_log_moles.rs → build_solve_contract()]
    ▼
┌─────────────────────────────────────────────────────┐
│  SolverPolicy::Cascade                              │
│  [equilibrium_solver_policy.rs]                     │
│                                                     │
│  for each backend in ordered_backends():             │
│      attempt = backend.solve(request)               │
│      [equilibrium_backend_adapter.rs]               │
│      if attempt успешен:                            │
│          validate_equilibrium_candidate()            │
│          [equilibrium_validation.rs]                │
│          if принят: → EquilibriumSolveCandidate     │
│      else: → SolverAttemptReport::Failed/Skipped    │
│                                                     │
│  Лучший кандидат → EquilibriumSolution              │
└─────────────────────────────────────────────────────┘
    │
    ▼
EquilibriumSolution
    │  log_moles, moles, validation report
    │  [equilibrium_problem.rs → EquilibriumSolution]
    ▼
MultiphaseEquilibriumSolution
    │  metadata, phase_totals, mole_fractions,
    │  phase_statuses, solve_report, keq_validation
    │  [phase_equilibrium_solution.rs]
```

### 3.2 Путь для однофазной газовой системы (исторический)

```
SubsData (список веществ)
    │  [User_substances.rs]
    ▼
gas_solver() / gas_solver_from_elements()
    │  search_substances() → parse_all_thermal_coeffs() →
    │  extract_all_thermal_coeffs(T) → calculate_elem_composition()
    │  → collect_gibbs_functions()
    │  [equilibrium_workflows.rs]
    ▼
EquilibriumLogMoles
    │  set_problem(n0, phases, P)
    │  [equilibrium_log_moles.rs]
    ▼
solve()
    │  → solve_candidate_from_seed()
    │  → каскад бэкендов
    │  → validate_equilibrium_candidate()
    │  [equilibrium_log_moles.rs, equilibrium_solver_policy.rs,
    │   equilibrium_backend_adapter.rs, equilibrium_validation.rs]
    ▼
EquilibriumSolution (внутренние поля solution, moles)
```

### 3.3 Путь с управлением фазами (Phase Control)

```
EquilibriumLogMoles::solve_with_phase_control()
    │  [equilibrium_workflows.rs]
    │
    │  1. Начальный PhaseSet из политики
    │  2. Цикл (max_phase_iterations):
    │     │
    │     ├── solve_fixed_active_set_candidate(active_mask)
    │     │   └── ActiveSetProjection::build()
    │     │       [equilibrium_active_set.rs]
    │     │       └── project_log_moles() → локальный решатель
    │     │       └── validate_element_totals_representable()
    │     │       └── solve_candidate_from_seed()
    │     │       └── scatter_log_moles() → глобальные координаты
    │     │
    │     ├── compute_phase_totals()
    │     │   [equilibrium_workflows.rs]
    │     ├── compute_phase_stability_reports()
    │     │   └── элементные потенциалы через SVD
    │     │   [equilibrium_workflows.rs]
    │     │
    │     ├── PhaseManager::classify_phases_at_temperature()
    │     │   └── PhaseTransitionPlan:
    │     │       ├── Deactivate → deactivate_phases_seed_only()
    │     │       ├── Activate → seed_activated_phase()
    │     │       ├── Hold → публикация с предупреждением
    │     │       └── NoTransition → публикация, выход
    │     │   [equilibrium_workflows.rs → PhaseManager]
    │     │
    │     └── reject_repeated_phase_set() (детекция циклов)
    │         [equilibrium_workflows.rs]
    │
    └── PhaseControlDidNotConverge (если лимит итераций)
```

### 3.4 Путь для температурной серии

```
solve_for_T_range(T_start, T_end, T_step)
    │  [equilibrium_log_moles.rs]
    │
    ├── Последовательный режим:
    │   for T in range:
    │       обновить gibbs(T)
    │       solve() с continuation seed (предыдущее решение)
    │       сохранить TemperatureSolveSnapshot
    │
    ├── Параллельный режим (par2):
    │   TemperatureWorkerSeed → rayon::par_iter()
    │   каждая точка независимо
    │
    └── postprocess_temperature_series()
        └── PCHIP интерполяция, передискретизация
        [equilibrium_temperature_postprocessing.rs]
```

### 3.5 Путь независимой K_eq валидации

```
EquilibriumSolution (каноническое)
    │  [equilibrium_problem.rs]
    ▼
EquilibriumConstantProblem
    │  (одна реакция, закон действующих масс)
    │  [equilibrium_constant_problem.rs]
    ▼
EquilibriumConstantSolver::solve()
    │  safeguarded Newton в пространстве extent
    │  bracket search, feasibility bounds
    │  [equilibrium_constant_solver.rs]
    ▼
EquilibriumConstantSolveResult
    │  [equilibrium_constant_solver.rs]
    ▼
EquilibriumConstantCrossValidationStatus
    │  сравнение: species_mole_delta, fraction_delta, total_Gibbs_delta
    │  [equilibrium_constant_cross_validation.rs]
```

---

## 4. Архитектурно-математические особенности

### 4.1 Log-mole формулировка

**Где:** [`equilibrium_log_moles.rs`](equilibrium_log_moles.rs), функции `evaluate_equilibrium_logmole_residual`, `evaluate_equilibrium_logmole_jacobian`

**Зачем:** Работа с малыми концентрациями (trace species) в обычных молях приводит к плохой обусловленности. Переход к `y_i = ln(n_i)` линеаризует область допустимых значений: вместо `n_i ≥ 0` получаем `y_i ∈ ℝ`.

**Как:** Невязка строится в пространстве реакционных координат (reaction extents). Система квадратная: `r` уравнений сродства + `e` уравнений элементного баланса, где `r = m - e` (m — число веществ, e — число элементов).

```math
f_k(y) = Σ_i ν_{ik} · y_i - Σ_j Δn_{kj} · ln(N_j) + Σ_j Δn_{kj} · φ_j + ln(K_k) = 0
f_{r+el}(y) = Σ_i a_{i,el} · exp(y_i) - b_el = 0
```

### 4.2 Символьный и численный стек

**Где:**
- Символьный: [`equilibrium_rst_backend.rs`](equilibrium_rst_backend.rs), [`equilibrium_workflows.rs`](equilibrium_workflows.rs) → `multiphase_equilibrium_residual_generator_sym`
- Численный: [`equilibrium_nonlinear.rs`](equilibrium_nonlinear.rs) → `LMSolver`, `NRSolver`, `TrustRegionSolver`
- Legacy адаптер: [`equilibrium_legacy_backend.rs`](equilibrium_legacy_backend.rs)
- Единый интерфейс: [`equilibrium_backend_adapter.rs`](equilibrium_backend_adapter.rs)

**Зачем:** Два независимых стека реализуют одну и ту же математическую формулировку. Символьный стек (RustedSciThe) строит residual и Jacobian как символьные выражения, затем лямбдифицирует их в численные функции. Численный стек (Legacy) использует hand-written Jacobian.

**Как:**
1. **Символьный путь:** `RstSymbolicThermochemistry::from_solver()` → `prepare_rst_symbolic_problem()` → `SymbolicNonlinearProblem` → RST решатели (6 методов: LM, MinpackLM, NielsenLM, TrustRegionLM, PowellDogleg, DampedNewton)
2. **Численный путь:** `solve_legacy_backend()` → `LMSolver`/`NRSolver`/`TrustRegionSolver` с hand-written residual и Jacobian

Выбор между стеками — через `SolverPolicy::RustedSciThe` vs `SolverPolicy::Legacy`.

### 4.3 Масштабирование (Row scaling и Variable scaling)

**Где:** [`equilibrium_problem.rs`](equilibrium_problem.rs) → `ResidualScalingContract`, `VariableScalingContract`

**Зачем:** Уравнения элементного баланса и реакции имеют разные порядки величин. Без масштабирования численный решатель может «видеть» только элементный баланс и игнорировать сродство реакций.

**Row scaling** — масштабирование строк residual и Jacobian:
```rust
ResidualScalingContract {
    scale: Vec<f64>,  // положительные множители для каждой строки
}
```
Применяется через `scale_residual_rows()` и `scale_jacobian_rows()`.

**Variable scaling** — масштабирование координат итерации:
```rust
VariableScalingContract {
    scale: Vec<f64>,  // положительные множители для каждой переменной
}
```
Применяется через `apply_iterate()` (деление на scale) и `unscale_iterate()` (умножение).

Два типа масштабирования разделены на уровне типов, чтобы их нельзя было перепутать.

### 4.4 Предотвращение нефизических значений

**Где:** Рассредоточено по всему модулю.

**Механизмы:**

1. **Trace floor** — `DEFAULT_TRACE_MOLE_FLOOR = 1e-30`. Нулевые начальные моли заменяются на положительный trace floor перед взятием логарифма. Политика: `TraceSpeciesSeedPolicy::Absolute` или `RelativeToLargestInitialMole`.

2. **LogMolesInitialGuess** — проверяет конечность и неотрицательность входных молей. Отрицательные значения отвергаются, а не обрезаются.

3. **EquilibriumConditions::new()** — проверяет, что T, P, P0 конечны и строго положительны.

4. **EquilibriumSolution::new()** — проверяет, что все log-moles конечны, все moles конечны и строго положительны.

5. **Feasibility check в решателях** — `max_step_moles_nonnegative()` ограничивает шаг NR, чтобы моли не стали отрицательными.

6. **Валидация кандидата** — `validate_equilibrium_candidate()` проверяет:
   - Конечность всех значений
   - Элементный баланс (независимо от residual)
   - Норму невязки
   - Минимальные моли

7. **Phase control** — `validate_phase_set_candidate()` проверяет, что неактивные фазы не содержат значимых молей.

### 4.5 SVD реакционный базис

**Где:** [`equilibrium_nonlinear.rs`](equilibrium_nonlinear.rs) → `compute_reaction_basis()`, [`equilibrium_reaction_basis.rs`](equilibrium_reaction_basis.rs) → `ValidatedReactionBasis`

**Зачем:** Из элементного состава системы нужно выделить независимые реакции. SVD разложение матрицы элементного состава даёт базис нуль-пространства, который и является набором независимых реакций.

**Как:**
1. SVD матрицы элементного состава A (m × e)
2. Ранг = число сингулярных чисел > tolerance
3. Число реакций r = m - rank
4. Базис нуль-пространства → матрица стехиометрии (m × r)
5. Gram-Schmidt достройка, если SVD не вернул полный базис

`ValidatedReactionBasis` дополнительно проверяет элементную сохранность: `A^T · N ≈ 0`.

### 4.6 Управление фазами: гистерезис, No flip-flop, bounded outer loop

**Где:** [`equilibrium_workflows.rs`](equilibrium_workflows.rs) → `PhaseManager`, `solve_with_phase_control()`

**Зачем:** В многофазной системе нужно решить, какие фазы активны. Простое правило «активна, если n > 0» приводит к флип-флопу (фаза появляется-исчезает на соседних итерациях).

**Гистерезис** — две границы: `dg_create` (создание) и `dg_keep` (удержание). Фаза создаётся только если `driving_force < dg_create`, удаляется только если `driving_force > dg_keep`. Между ними — гистерезисная петля.

```rust
enum PhaseHysteresisPolicy {
    TemperatureScaled { create_rt_factor, keep_rt_factor },
    Explicit { dg_create, dg_keep },
}
```

По умолчанию: `create_rt_factor = -1e-6`, `keep_rt_factor = 1e-8`. При T=1000K: `dg_create = -0.008314`, `dg_keep = 0.00008314`.

**No flip-flop:**
- Один переход за итерацию outer loop (деактивация имеет приоритет)
- `reject_repeated_phase_set()` — детекция циклов через HashSet посещённых PhaseSet
- `max_phase_iterations = 16` — bounded outer loop

**Bounded outer loop:**
```
for iteration in 0..max_phase_iterations:
    solve_fixed_active_set()
    compute_phase_stability()
    classify_phases() → PhaseTransitionPlan
    match:
        Deactivate → seed trace, continue
        Activate → seed 1e-8 от system total, continue
        Hold → публикация с Hold статусом
        NoTransition → публикация, return
return PhaseControlDidNotConverge
```

### 4.7 Phase stability (проверка стабильности фаз)

**Где:** [`equilibrium_workflows.rs`](equilibrium_workflows.rs) → `compute_phase_stability_reports()`

**Зачем:** Определить, должна ли неактивная фаза появиться, или активная — исчезнуть.

**Как:**
1. Для чистой конденсированной фазы: вычисляются элементные потенциалы λ через SVD активных веществ: `A_active^T · λ = μ_active`
2. Driving force: `ΔG = μ_candidate - Σ_j a_{candidate,j} · λ_j`
3. Если `ΔG < 0` — фаза термодинамически выгодна (должна появиться)
4. Если `ΔG > 0` — фаза нестабильна (должна исчезнуть)

Для газовой фазы используется модель `FixedIdealGas` (опорные потенциалы, но фиксированная активность). Многокомпонентные растворы пока не поддерживаются (требуют tangent-plane minimization).

### 4.8 Element-validation after activation

**Где:** [`equilibrium_active_set.rs`](equilibrium_active_set.rs) → `validate_element_totals_representable()`

**Зачем:** После активации/деактивации фазы нужно убедиться, что новый набор активных веществ может представить исходный элементный инвентарь системы.

**Как:**
1. Берётся транспонированная матрица элементного состава активных веществ: `A_active^T`
2. Решается СЛАУ: `A_active^T · x = b` (b — полные элементные моли)
3. Вычисляется невязка: `||A_active^T · x - b||`
4. Если невязка > tolerance — активный набор не может представить элементный состав, фазовая конфигурация отвергается

Это предотвращает ситуацию, когда деактивация фазы удаляет единственное вещество, содержащее некоторый элемент.

### 4.9 Независимая K_eq валидация

**Где:**
- [`equilibrium_constant_problem.rs`](equilibrium_constant_problem.rs) — постановка задачи
- [`equilibrium_constant_solver.rs`](equilibrium_constant_solver.rs) — решатель
- [`equilibrium_constant_validation.rs`](equilibrium_constant_validation.rs) — отчёт
- [`equilibrium_constant_cross_validation.rs`](equilibrium_constant_cross_validation.rs) — сравнение

**Зачем:** Основной решатель использует residual/Jacobian формулировку. Независимый K_eq решатель использует закон действующих масс — полностью другую математическую формулировку. Если два независимых решателя дают одинаковый ответ, доверие к решению выше.

**Как:**
1. Строится `EquilibriumConstantProblem` для одной реакции
2. Решается safeguarded Newton в пространстве extent ξ:
   - `f(ξ) = ln(Q(ξ)) - ln(K(T))`
   - Поиск скобки (bracket search) с учётом границ осуществимости
   - Safeguarded Newton с finite-difference производной
3. Результат: `EquilibriumConstantSolveResult` (extent, moles, validation)
4. Cross-validation: сравнение с каноническим решением по:
   - `max_abs_species_mole_delta`
   - `max_abs_species_fraction_delta`
   - `max_abs_total_gibbs_delta`

Режимы: `Off`, `WhenApplicable` (только для однореакционных систем), `Required`.

### 4.10 Транзакционность публикации состояния

**Где:** [`equilibrium_log_moles.rs`](equilibrium_log_moles.rs) → `publish_solve_candidate()`, `EquilibriumSolveCandidate`

**Зачем:** В phase-control outer loop промежуточные решения не должны перезаписывать последнее опубликованное состояние. Только когда outer loop сошёлся, решение публикуется.

**Как:**
1. `solve_fixed_active_set_candidate()` возвращает `EquilibriumSolveCandidate` — временный результат
2. Phase control анализирует кандидата, решает, нужен ли переход
3. Только при `NoTransition` или `Hold` вызывается `publish_solve_candidate()`, которая:
   - Записывает `self.solution`, `self.moles`
   - Сохраняет `last_validation_report`, `last_solve_report`
   - Сохраняет `last_keq_validation_status`
   - Сохраняет `last_phase_control_report`

### 4.11 Политика решателя (SolverPolicy) и каскад

**Где:** [`equilibrium_solver_policy.rs`](equilibrium_solver_policy.rs), [`equilibrium_backend_adapter.rs`](equilibrium_backend_adapter.rs)

**Зачем:** Разные численные методы имеют разные профили сходимости. Каскад пробует методы по порядку, пока один не сойдётся.

**Как:**
1. `SolverPolicy::Single(backend)` — ровно одна попытка
2. `SolverPolicy::Cascade(vec)` — упорядоченный список, дубликаты удаляются
3. `SolverCascadeBudget` — ограничения:
   - `max_attempts`: максимум запущенных бэкендов
   - `max_iterations_per_attempt`: итераций на один бэкенд
   - `max_total_iterations`: суммарно на весь каскад
4. Каждый бэкенд возвращает `SolverAttemptReport` с:
   - `outcome`: Accepted / Failed / RejectedCandidate / Skipped
   - `failure_kind`: Solver / Backend / ResidualEvaluation / JacobianEvaluation
   - `metrics`: termination, iterations, residual_evaluations, elapsed_millis

Рекомендуемый RST каскад: LM → MinpackLM → NielsenLM → TrustRegionLM → PowellDogleg → DampedNewton.

### 4.12 Температурные серии и постобработка

**Где:** [`equilibrium_temperature_postprocessing.rs`](equilibrium_temperature_postprocessing.rs)

**Зачем:** После расчёта равновесия в диапазоне температур нужно построить гладкие кривые для визуализации.

**Как:**
1. `TemperatureSweepSeries` — сырые данные (температура, метки, строки значений)
2. `TemperaturePostprocessingPolicy`:
   - `grid`: RawOnly / Uniform { points } / Explicit(vec)
   - `interpolation`: Linear / Log (для положительных величин)
3. `postprocess_temperature_series()` — PCHIP интерполяция, гарантирующая монотонность

---

## 5. Навигационная карта

### 5.1 С чего начать чтение

**Если вы хотите понять, как работает расчёт равновесия:**

1. **Начните с [`equilibrium_problem.rs`](equilibrium_problem.rs)** — это входная граница. Поймите `EquilibriumProblem`, `EquilibriumConditions`, `EquilibriumSolution`, `LogMolesInitialGuess`.

2. **Затем [`equilibrium_log_moles.rs`](equilibrium_log_moles.rs)** — основной оркестратор. Сосредоточьтесь на:
   - `EquilibriumLogMoles::solve()` — однотемпературный расчёт
   - `EquilibriumLogMoles::solve_candidate_from_seed()` — ядро каскада
   - `evaluate_equilibrium_logmole_residual()` — математическая формулировка

3. **Потом [`equilibrium_workflows.rs`](equilibrium_workflows.rs)** — управление фазами:
   - `PhaseManager`, `PhaseSet`, `PhaseTransitionPlan`
   - `solve_with_phase_control()` — bounded outer loop
   - `compute_phase_stability_reports()` — проверка стабильности

**Если вы хотите добавить новый бэкенд:**

1. [`equilibrium_backend_adapter.rs`](equilibrium_backend_adapter.rs) — trait `EquilibriumNonlinearBackend`
2. [`equilibrium_solver_policy.rs`](equilibrium_solver_policy.rs) — добавить вариант в `SolverBackend`
3. [`equilibrium_rst_backend.rs`](equilibrium_rst_backend.rs) или новый файл адаптера

**Если вы хотите понять валидацию:**

1. [`equilibrium_validation.rs`](equilibrium_validation.rs) — каноническая проверка
2. [`equilibrium_constant_solver.rs`](equilibrium_constant_solver.rs) — независимый K_eq решатель
3. [`equilibrium_constant_cross_validation.rs`](equilibrium_constant_cross_validation.rs) — сравнение

### 5.2 Карта файлов

| Файл | Назначение | Ключевые типы |
|------|------------|---------------|
| [`equilibrium_ids.rs`](equilibrium_ids.rs) | Типобезопасные идентификаторы | `SpeciesId`, `ElementId`, `PhaseIndex`, `ReactionId` |
| [`equilibrium_component.rs`](equilibrium_component.rs) | Фазово-квалифицированный компонент | `EquilibriumComponentDescriptor` |
| [`equilibrium_problem.rs`](equilibrium_problem.rs) | Граница задачи | `EquilibriumProblem`, `PreparedEquilibriumProblem`, `EquilibriumSolution`, `LogMolesInitialGuess`, `ResidualScalingContract`, `VariableScalingContract` |
| [`equilibrium_log_moles.rs`](equilibrium_log_moles.rs) | Основной оркестратор (3881 строка) | `EquilibriumLogMoles`, `EquilibriumSolverSettings`, `Solvers`, `SolverParams`, `EquilibriumSolveCandidate`, `TemperatureWorkerSeed` |
| [`equilibrium_nonlinear.rs`](equilibrium_nonlinear.rs) | Численные решатели | `LMSolver`, `NRSolver`, `TrustRegionSolver`, `ReactionBasis`, `ReactionExtentError` |
| [`equilibrium_workflows.rs`](equilibrium_workflows.rs) | Управление фазами, stability, convenience (2043 строки) | `PhaseManager`, `PhaseSet`, `PhaseTransitionPlan`, `PhaseStabilityReport`, `PhaseControlledSolveReport`, `MultiphaseAcceptanceReport`, `gas_solver()` |
| [`equilibrium_active_set.rs`](equilibrium_active_set.rs) | Проекция активного набора фаз | `ActiveSetProjection` |
| [`equilibrium_activity.rs`](equilibrium_activity.rs) | Модели активности фаз | `PhaseActivityModel` |
| [`equilibrium_validation.rs`](equilibrium_validation.rs) | Каноническая проверка кандидатов | `EquilibriumCandidateReport`, `EquilibriumAcceptanceCriteria` |
| [`equilibrium_solver_policy.rs`](equilibrium_solver_policy.rs) | Политика выбора бэкенда | `SolverPolicy`, `SolverBackend`, `SolverCascadeBudget`, `SolverAttemptReport` |
| [`equilibrium_backend_adapter.rs`](equilibrium_backend_adapter.rs) | Единый интерфейс бэкенда | `EquilibriumNonlinearBackend`, `BackendSolveRequest` |
| [`equilibrium_rst_backend.rs`](equilibrium_rst_backend.rs) | Адаптер RustedSciThe (символьный) | `RstPreparedProblem`, `RustedSciTheSolver`, `RustedSciTheSolveOutcome` |
| [`equilibrium_legacy_backend.rs`](equilibrium_legacy_backend.rs) | Адаптер Legacy (численный) | `solve_legacy_backend()` |
| [`equilibrium_constant_problem.rs`](equilibrium_constant_problem.rs) | K_eq постановка задачи | `EquilibriumConstantProblem` |
| [`equilibrium_constant_solver.rs`](equilibrium_constant_solver.rs) | K_eq решатель | `EquilibriumConstantSolver`, `EquilibriumConstantSolveResult` |
| [`equilibrium_constant_validation.rs`](equilibrium_constant_validation.rs) | K_eq отчёт валидации | `EquilibriumConstantValidationReport` |
| [`equilibrium_constant_cross_validation.rs`](equilibrium_constant_cross_validation.rs) | Сравнение канонического и K_eq | `EquilibriumConstantCrossValidationStatus` |
| [`equilibrium_reaction_basis.rs`](equilibrium_reaction_basis.rs) | Типобезопасный реакционный базис | `ValidatedReactionBasis` |
| [`equilibrium_multiphase_domain.rs`](equilibrium_multiphase_domain.rs) | Многофазная доменная граница | `MultiphaseEquilibriumLayout`, `MultiphaseInitialComposition` |
| [`phase_equilibrium_problem.rs`](phase_equilibrium_problem.rs) | Мост к фазовой подсистеме | `PhaseEquilibriumMetadata`, `PhaseEquilibriumBuildRequest`, `PhaseEquilibriumSolutionBundle` |
| [`phase_equilibrium_solution.rs`](phase_equilibrium_solution.rs) | Иммутабельный многофазный результат | `MultiphaseEquilibriumSolution` |
| [`phase_equilibrium_workflow.rs`](phase_equilibrium_workflow.rs) | Публичный фасад fixed-P,T | `solve_resolved_pt()` |
| [`equilibrium_temperature_postprocessing.rs`](equilibrium_temperature_postprocessing.rs) | Постобработка температурных серий | `TemperaturePostprocessingPolicy`, `TemperatureSweepSeries` |
| [`easy_equilibrium.rs`](easy_equilibrium.rs) | Упрощённый интерфейс (одна реакция) | `EasyEquilibrium` |
| [`NR_Legacy.rs`](NR_Legacy.rs) | Legacy Newton-Raphson | `NRSolver` (legacy) |

### 5.3 Типичные сценарии использования

**Сценарий 1: Однофазный газовый расчёт**
```rust
let mut solver = gas_solver(
    vec!["CO2".into(), "CO".into(), "O2".into()],
    1000.0, 101325.0, Solvers::LM, Some("info"), true
)?;
solver.set_problem(vec![1.0, 0.0, 0.0], phases, 101325.0)?;
solver.solve()?;
```

**Сценарий 2: Многофазный расчёт через новый API**
```rust
let request = ResolvedPhaseEquilibriumRequest::new(
    &resolved_system,
    EquilibriumConditions::new(1000.0, 101325.0, 101325.0)?,
    initial_composition,
);
let solution = solve_resolved_pt(request)?;
```

**Сценарий 3: Расчёт с управлением фазами**
```rust
solver.phase_manager = PhaseManager::new(1e-30, -1e-6, 1e-8);
solver.solve_with_phase_control()?;
```

**Сценарий 4: Температурная серия**
```rust
solver.solve_for_T_range(500.0, 2500.0, 100.0)?;
let processed = postprocess_temperature_series(
    &solver.temperature_solutions,
    TemperaturePostprocessingPolicy::default(),
)?;
```

---

## Приложение: Диаграмма зависимостей модулей

```
equilibrium_ids.rs
    │
    ├── equilibrium_component.rs
    │       │
    │       └── equilibrium_problem.rs
    │               │
    │               ├── equilibrium_log_moles.rs
    │               │       │
    │               │       ├── equilibrium_nonlinear.rs
    │               │       │       └── (LMSolver, NRSolver, TrustRegionSolver)
    │               │       │
    │               │       ├── equilibrium_workflows.rs
    │               │       │       ├── equilibrium_active_set.rs
    │               │       │       ├── equilibrium_activity.rs
    │               │       │       └── (PhaseManager, PhaseSet, stability)
    │               │       │
    │               │       ├── equilibrium_backend_adapter.rs
    │               │       │       ├── equilibrium_rst_backend.rs
    │               │       │       └── equilibrium_legacy_backend.rs
    │               │       │
    │               │       ├── equilibrium_solver_policy.rs
    │               │       ├── equilibrium_validation.rs
    │               │       ├── equilibrium_constant_solver.rs
    │               │       │       ├── equilibrium_constant_problem.rs
    │               │       │       └── equilibrium_constant_validation.rs
    │               │       ├── equilibrium_constant_cross_validation.rs
    │               │       └── equilibrium_temperature_postprocessing.rs
    │               │
    │               └── phase_equilibrium_problem.rs
    │                       │
    │                       ├── equilibrium_multiphase_domain.rs
    │                       └── phase_equilibrium_solution.rs
    │
    └── phase_equilibrium_workflow.rs
            │
            └── (solve_resolved_pt — публичный вход)