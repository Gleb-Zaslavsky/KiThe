# Chemical equilibrium calculator

This form creates one canonical `EquilibriumCalculator` request. It does not
contain a separate GUI solver. Results are immutable accepted snapshots.

## Setup

`P,T = const` solves composition at fixed pressure and temperature. `P,H =
const` solves composition and temperature at fixed pressure and total
enthalpy. `Point` solves one state. `Range` is available for P,T and keeps the
entered direction, reusing formulation data and the last accepted state.

`Pressure [Pa]` is the system pressure and must be finite and positive.
`Reference pressure [Pa]` is the gas standard-state pressure in the ideal-gas
activity law; it is not automatically the system pressure and must match the
library convention. `Temperature [K]` is the P,T temperature.

For a P,T range, `Start`, `End`, and `Points` define the grid; at least two
points are required. For P,H, `Target total enthalpy [J]` is extensive enthalpy
of the complete inventory, not J/mol. `Lower bound`, `Upper bound`, and
`Initial seed [K]` define and seed the scalar temperature search. The seed must
be inside the bracket, and the bracket must stay inside available data.

`Explicit species` accepts exact records. `Search by elements` discovers
candidates, but candidates must still be assigned to a phase. `Elements` is the
requested element set; `Exact` requires the same set while `Subset` allows a
candidate subset. State filters and candidate limits make catalog searches
deterministic.

Both input modes are prepared through the public typed
`ChemEquilibrium::prelude` facade. Explicit species provide the molecular
initial state; element mode provides the conserved elemental inventory, while
candidate search only selects the real phase-qualified species universe.

For every phase, `Phase` is its stable identity, `Physical state` is gas/liquid/
solid/condensed semantics, and `Model` selects ideal gas, ideal
solution, or pure condensed behavior. Incompatible state/model combinations are
rejected. `Substance` is a library identifier and `Initial moles` is a
phase-qualified amount. `Add phase`, `Add component`, and `Remove` edit only
the document and require preparation again.

`Validate document` checks field-local structure. `Prepare canonical request`
maps validated state to the facade without running a solve. `Run prepared
request` starts the worker. `Cancel` prevents late or stale output from
replacing accepted state.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Problem, point/range, pressure and temperature | Positive Pa and K; `Point` is the default | Always; range is P,T only | Request-invalidating |
| P,H target and temperature bracket | Extensive J; finite K bracket and seed inside it | P,H only | Request-invalidating |
| Species, elements, phases and initial moles | Exact identifiers; finite non-negative mol | According to inventory mode | Request-invalidating |
| Validate / Prepare / Run / Cancel | Validate first; preparation creates an immutable request | Action availability follows document and worker state | Preparation-only / runtime |

## Phase control

`Fixed declared phases` keeps exactly the declared phase set. `Bounded phase
control` enables the outer lifecycle, which may create or remove candidates
after accepted TPD/stability checks.

`Phase epsilon [mol]` is the lifecycle presence threshold. `Creation driving
force` controls appearance of an absent phase; `Keep driving force` controls
retention of an active phase. Their separation is hysteresis and prevents
create/remove chatter. `Maximum phase iterations` bounds lifecycle work; an
unverified intermediate state is never published.

`Positive inventory` starts from phases with positive amount. `All declared
candidates` also tests zero-inventory declared phases, which is useful for
appearance tests. These are lifecycle policies, not new physical models.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Fixed / bounded lifecycle | Fixed declared phases is the default | Always | Request-invalidating |
| Phase epsilon | Positive mol threshold | Bounded lifecycle | Request-invalidating |
| Create / keep driving force | Finite J/mol thresholds; separated values form hysteresis | Bounded lifecycle | Request-invalidating |
| Initial phase set and iteration budget | Positive inventory; finite positive iteration count | Bounded lifecycle | Request-invalidating |

## Libraries

`Engine default` delegates priority to the shared repository. `Explicit policy`
uses the ordered `Priority libraries` list and closed `Permitted libraries` set.
`Offline mode` prevents network access. `Allow online NIST fallback` is an
explicit opt-in and never silently replaces a local record.

`Library catalog status` reports loading, loaded, unavailable, or failed state
of the shared catalog. `Lookup instruction` describes the requested source
before resolution. Actual record identity, library, version, and provenance
are authoritative only in accepted Results.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Lookup policy | Engine default unless explicit policy is selected | Always | Request-invalidating |
| Priority and permitted libraries | Ordered priority and closed allow-list | Explicit policy | Request-invalidating |
| Offline / NIST fallback | Offline is safe default; NIST requires opt-in | Always | Request-invalidating |
| Load catalog | Shared local repository action | Catalog is not ready | Preparation-only |

## Numerics

`Production default cascade` uses the configured backend cascade. `Single backend` is
useful for diagnosis and has no fallback. `Custom cascade` is for controlled
comparisons and must contain unique supported backends.

`P,H route` offers `Auto`, `Monolithic`, and `Nested temperature` where
supported. `Tolerance` controls numerical acceptance; `Maximum iterations`
limits work. Empty overrides use production policy. `Enable scaling` changes
solver coordinates, not physical units. `Trace seed policy` sets positive
log-mole initial values and is not extra inventory. `Cascade total iterations`
limits total fallback work while retaining attempt evidence.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Solver policy | Production cascade by default | Always | Request-invalidating |
| Backend and custom cascade | Supported unique backends | Single backend / custom cascade | Request-invalidating |
| P,H route | Auto by default | P,H only | Request-invalidating |
| Tolerance, iteration and scaling overrides | Empty uses production defaults; finite positive overrides | Advanced numerics | Request-invalidating |
| Trace-seed policy | No physical inventory is added | Advanced numerics | Request-invalidating |

## Output

`Collect timing` records stage and point timings. `Phase lifecycle trace` keeps
TPD, stability, activation, deactivation, and transition decisions for bounded
control. `Maximum retained lifecycle events` bounds trace memory while keeping
recent events. `Range lifecycle trace` selects endpoints, transitions only, or
every point.

`K_eq validation` enables an independent small-system validation route; it is
not a replacement for the general production solver. `Result basis` selects
moles, mole fractions, or phase totals for derived display series.

`Result table` chooses detailed moles plus fractions or compact moles only.
`Hide below mole fraction` filters rows in the visible table; zero shows all
components. Both are display-only and do not change the snapshot.

`Plot target` selects embedded plot, KiThePlot, both, or none. `Y scale` selects
linear or log10 presentation. `PCHIP display resampling` is available for P,T
ranges; its point count, interpolation space, and clamp affect display data
only, never solver points.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Timing and lifecycle trace | Disabled by default; retained-event cap bounds memory | Diagnostics options | Runtime diagnostics |
| K_eq validation | Disabled by default; small applicable systems only | Diagnostics options | Request-invalidating |
| Table, basis, threshold and visible series | Presentation defaults; finite fraction threshold | Accepted results | Display-only |
| Plot target and scale | No plot / linear are safe defaults | Accepted results | Display-only |
| PCHIP resampling | Disabled; display-point count applies to P,T range | P,T range results | Display-only |

## Results

Each accepted point is an immutable accepted snapshot and is read-only. `Solve diagnostics` contains residuals,
element conservation, backend attempts, acceptance, fallback reasons, timing,
and validation evidence. `Lookup provenance` shows the resolved source for each
component. `Phase lifecycle` shows TPD and transition evidence when retained.
Range summaries show direction, accepted points, formulation builds/reuses,
normalization recoveries, transitions, and timing. Failed points are not shown
as accepted points, and plotting never runs the solver again.

### Control contract

| Control | Domain / effective default | Visible when | Change class |
| --- | --- | --- | --- |
| Accepted snapshot | Read-only accepted physical state | A solve has succeeded | Display-only |
| Diagnostics, provenance and lifecycle evidence | Read-only evidence retained by the request | Matching collection was enabled | Display-only |
| Range selection and plots | Derived from accepted range points only | Accepted P,T range | Display-only |
