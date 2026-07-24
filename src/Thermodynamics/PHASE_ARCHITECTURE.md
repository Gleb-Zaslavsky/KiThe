# Phase Subsystem Architecture

This document describes the \`phase_*\` subsystem in the thermodynamics layer.
It covers the boundary between thermodynamic-data lookup (\`SubsData\`) and
higher-level equilibrium or solver code. It does not define an equilibrium
algorithm or silently select a physical activity model.

## Core Terms

| Term | Meaning in this codebase |
| --- | --- |
| **Transactional operation** | A multi-phase operation works on cloned payloads and staged results, then commits only if every phase and validation step succeeds. |
| **Cache** | Rebuildable derived state: numeric Gibbs/entropy values, closures, symbolic expressions, and evaluation context. It is not source library data. |
| **Snapshot** | A typed, revision-aligned read-only cache view with the layout and conditions that make results valid. |
| **Cache bundle** | \`PhaseThermoCacheBundle\`: one phase's cache grouped by property family and representation. |
| **Layout** | \`SystemLayout\`, the sole canonical ordering of phase-qualified solver components. |
| **Resolution** | Fallible mapping from declarative \`PhaseSpec\` to selected library records, layout, and provenance. |

## Purpose and Boundary

The subsystem sits between two layers:

\`\`\`text
Thermodynamic libraries / SubsData
    -> lookup, parsing, coefficients, transport and thermo records

phase_* subsystem
    -> typed phase declarations, deterministic component order,
       pure evaluation, transactional transitions, cache validity

Chemical-equilibrium and reactor consumers
    -> equations, equilibrium constraints, numerical solvers
\`\`\`

Its purpose is to give upper layers coherent phase data: which components
exist, which solver order they use, which thermodynamic functions belong to
them, and under which conditions a result was obtained. It does not replace a
physical model or an equilibrium solver.

## Structures and Dependencies

\`\`\`text
PhaseSpec --------------------------> ResolvedPhaseSystem
  |                                      |
  | PhaseId, components,                 | resolved SubsData payloads
  | physical state, model                | SystemLayout + provenance
  v                                      v
phase_domain.rs                     PhaseSystem
                                         |
              +--------------------------+--------------------------+
              |                          |                          |
              v                          v                          v
        phase_layout.rs            phase_evaluation.rs       phase_cache.rs
        canonical order            typed pure requests       bundles and snapshots
              |                          |                          |
              +--------------------------+--------------------------+
                                         |
                                         v
                    phase_operations.rs + phase_property_builders.rs
                    staged fallible work, no state publication
\`\`\`

\`PhaseSystem\` is the only orchestration module. It owns phase \`SubsData\`,
cache bundles, revisions, resolved layout, and numeric-evaluation context.
Extracted modules cannot publish state independently: they either describe
data, build a temporary result, or expose a read-only view.

This protects against partial publication. If one phase has obtained new
coefficients or closures and processing of the next phase fails, external code
must not observe a mixed state.

## Data Lifecycle

\`\`\`text
1. Declare: PhaseSpec + lookup policy
        |
2. Resolve: ThermoRepository / SubsData lookup
        |
        +--> ResolvedPhaseSystem + provenance + SystemLayout
        |
3. Install: PhaseSystem installs payloads and invalidates derived state
        |
4. Evaluate/build: numeric request OR symbolic/closure staged construction
        |
5. Validate: phases, component order, finite values, conditions, revision
        |
6. Commit: one PhaseSystem publication and revision update
        |
7. Read: borrowed view or typed snapshot aligned to current layout
\`\`\`

Revision and cache context are part of the contract, rather than debugging
metadata. After a payload, layout, or symbolic expression changes, an old
result must not appear current.

## Why Components Are Phase-Qualified

\`PhaseComponentId\` includes both phase and substance. \`H2O(g)\` and \`H2O(l)\`
may have the same elemental composition, but they are different solver
unknowns: their standard states, chemical potentials, and mixing models can
differ. Storing them only by substance name would merge or overwrite a
component incorrectly.

\`SystemLayout\` defines the stable order of such components. \`HashMap\` remains
appropriate for lookup and compatibility with library data, but never defines
solver-vector order. This prevents non-deterministic vectors, coefficient
permutations, and matrix errors.

## Transactional Mutation

\`\`\`text
visible PhaseSystem state
        |
        v
clone every phase payload
        |
        v
run a fallible operation for every phase
        |
   error? ---- yes ---> discard the whole stage
        |
        no
        v
validate results, request shape, finite values, and layout
        |
        v
one commit + revision update + cache invalidation/publication
\`\`\`

The same rule applies to numeric cache publication: the complete request,
phase and component coverage, finite values, and conditions are checked before
a cache family and its context are replaced. An error therefore cannot leave a
new Gibbs cache with an old context, or the reverse.

## Cache, Bundle, and Snapshot

\`\`\`text
PhaseThermoCacheBundle for one phase
|
+-- Gibbs property cache
|   +-- numeric: substance -> f64
|   +-- functions: substance -> closure
|   \`-- symbolic: substance -> Expr
|
\`-- Entropy property cache
    +-- numeric: substance -> f64
    +-- functions: substance -> closure
    \`-- symbolic: substance -> Expr
\`\`\`

A number, a closure, and an \`Expr\` are not treated as independent global maps.
They are representations of one property in one phase. The bundle keeps this
relationship in the data model and reduces the risk of updating \`dG_sym\` while
forgetting to invalidate the corresponding numeric \`dG\`.

A snapshot adds validity conditions to cached data: revision, \`T\`, \`P\`, and
normalized composition for numeric evaluation. It can also project numeric
values into a dense vector in exact \`SystemLayout\` order for solver hot paths.

## Why Pure Evaluation Is Separate from Publication

\`phase_evaluation.rs\` defines and validates a typed request before a calculator
runs. \`phase_property_builders.rs\` creates symbolic expressions and closures
against cloned payloads. These steps can return \`SubsDataResult\`, but do not
mutate \`PhaseSystem\`.

This gives three properties:

1. A failure in one phase cannot corrupt another phase.
2. Computation can be tested without inspecting hidden cache state.
3. Exactly one location owns revision, invalidation, and publication policy.

This boundary is not abstraction for its own sake. It keeps mutable \`SubsData\`
where it is currently needed to acquire coefficient data without propagating
implicit mutations into upper layers.

## Public API and Legacy Boundary

External callers use \`OnePhase\`, \`PhaseOrSolution\`, typed specifications, and
read-only views. \`PhaseSystem\` is crate-private because direct mutation of its
maps would bypass invalidation and layout revisioning.

Legacy \`Fn(...) -> f64\` closures remain a compatibility representation, but
cannot return a typed error at invocation time. New numeric paths should prefer
the fallible typed request/evaluation API. Full migration of legacy equilibrium
consumers is intentionally deferred until work on \`ChemEquilibrium\`.

## Module Responsibilities

| Module | Responsibility | It must not own |
| --- | --- | --- |
| \`phase_domain.rs\` | Phase declarations, physical state/model, resolved records, provenance | cache mutation or evaluation policy |
| \`phase_layout.rs\` | Canonical ordered phase/component layout | database lifecycle or cache publication |
| \`phase_evaluation.rs\` | Typed numeric request validation and pure evaluation | cache ownership |
| \`phase_operations.rs\` | Transactional staging over cloned \`SubsData\` payloads | commit or revision policy |
| \`phase_property_builders.rs\` | Staged symbolic and closure construction | visible state publication |
| \`phase_cache.rs\` | Cache bundle types, views, snapshots, dense projections | mutable lifecycle ownership |
| \`phase_system.rs\` | Commit, invalidation, revisioning, cache lifecycle | a second public facade |
| facade modules | Ergonomic and compatibility-facing API | independent engine state |

## Rule for Future Changes

For a new physical model, property, or solver adapter, first decide whether it
is declarative domain data, pure evaluation, a staged transformation, or a
\`PhaseSystem\` lifecycle transition. Only the final category may mutate visible
phase state or publish caches. This keeps the subsystem deterministic,
testable, and resistant to partial updates.

