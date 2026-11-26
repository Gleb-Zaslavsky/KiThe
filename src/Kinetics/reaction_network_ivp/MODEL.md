# Reaction-Network IVP Model

The initial implementation solves a closed homogeneous constant-volume network
in physical time. The prescribed temperature is an input, not an integrated
state. There is no energy balance or thermochemical lookup.

## Coordinates And Units

- Independent coordinate `t`: seconds.
- State `C_i(t)`: species concentration in mol/m^3.
- SolidState state `alpha_j(t)`: dimensionless conversion, after the
  concentration coordinates in the prepared IVP.
- Prescribed `T(t)`: kelvin, positive wherever the model is evaluated.
- Forward reaction rate `r_j`: mol reaction/(m^3*s).
- Stoichiometric coefficient `nu[j,i]`: products minus reactants, dimensionless.
- State equation: `dC_i/dt = sum_j nu[j,i] * r_j(C,T,t)`.

Species and reactions retain input order and stable IDs. The internal symbols
`c0`, `c1`, ... follow species order, regardless of user-facing IDs. `t` is
the independent variable and `T` in a custom rate denotes the prescribed
temperature expression. The equations are prepared symbolically for RST; the
Jacobian is taken with respect to the full `[C, alpha]` state at fixed `t`.
Constant temperature is prepared with numeric parameter `p_temperature`.
Linear heating uses `T(t)=p_temperature+p_heating_rate*t`; `p_temperature`
is the value extrapolated to `t=0`, even if integration starts later.
`PreparedNetworkIvp::solve_with_temperature_binding` rebinds these numeric
values and optionally the initial concentrations without rebuilding the Expr
graph. A different temperature-law shape requires a new preparation.
`PreparedNetworkIvp::solve_with_rate_coefficients` rebinds the complete set of
numeric coefficients for `MassAction`, `EmpiricalMassAction`, and `SolidState`, keyed by
reaction ID. Constant coefficients bind one value; Arrhenius coefficients bind
the prefactor and activation energy together. The binding list may be supplied
in any order, but omissions, duplicates, unknown IDs, and changes of
coefficient kind are rejected before integration. A combined method can bind
temperature and rate coefficients in one fresh solve. Accepted bindings are
published only with an accepted result; failed requests do not change the
last-good result or accepted numeric values.
An empty species list is rejected. An empty reaction list is accepted as an
inert network, with one zero RHS per species. Species absent from all reactions
also have zero RHS. Repeated equations with distinct reaction IDs are distinct
stages and contribute independently.

## Implemented Rate Laws

For `MassAction`, the order of each reactant is its forward stoichiometric
coefficient:

`r_j = k_j(T) * product_reactants C_i^a[j,i]`.

`RateCoefficient::Constant(k)` takes the rate coefficient in units consistent
with that order. For `Arrhenius`,

`k(T) = A * exp(-E_a/(R*T))`, `R = 8.31446261815324 J/(mol*K)`.

`E_a` is J/mol. The units of `A` and `k` depend on total reaction order. A
negative finite activation energy is permitted as an empirical law. A zero
prefactor is permitted. Reverse reactions must be explicit independent stages.

`EmpiricalMassAction` uses explicit non-negative kinetic exponents instead of
reactant stoichiometric coefficients. These exponents alter only the rate;
the net species balance still uses products minus reactants. An empty order
list denotes a zero-order rate. Negative exponents can instead be expressed
as `CustomExpr`, with their singular-domain responsibility explicit.

`CustomExpr` is the *complete* forward rate. The calculator adds no hidden
mass-action or Arrhenius multiplier. It may use only `c0..cN`, `t`, and `T`.
Parsed strings are not part of this initial API; callers supply RST `Expr`.

`SolidState` reuses `KineticModel::get_expr()` inside the same vector IVP.
Each stage supplies capacity `q_j > 0` in mol reaction/m^3 and `alpha_j(0)`;
`alpha_j' = k_j(T) f_j(alpha_j)` and `r_j = q_j alpha_j'`. Neither capacity
nor conversion is inferred from concentration. Each net-consumed species
must initially cover the stage's remaining capacity:
`C_i(0) >= -nu[j,i] q_j (1-alpha_j(0))`. Future production is not credited.
Two stages cannot both net-consume the same species if one is SolidState.
Conversion stays in `[0,1]` and in the model's solver-admissible domain;
an endpoint is excluded if its state Jacobian is singular or its nonzero
rate would immediately drive conversion past one. Only F1/F2/F3 currently
admit both endpoints;
rates are nonnegative. There is no clipping, epsilon shift, or implicit
reverse conversion. Output carries conversions separately, in SolidState
reaction order. Rebinding concentrations retains `alpha_j(0)` and rechecks
capacity.

For `A -> 2B` and `B+C -> D`, the RHS is exactly

`A'=-r1`, `B'=2*r1-r2`, `C'=-r2`, `D'=r2`.

The current public facade exposes BDF, Radau, Backward Euler, and LSODE2 via
StiffAPI. Its default is BDF, ExprLegacy, Dense, sequential callbacks where
supported, with native telemetry and logging Off. Unsupported settings are
reported as typed KiThe or RST capability errors, not silently replaced.
`rtol` and `atol` control BDF, Radau, and LSODE2 integration. Backward Euler
uses its native Newton tolerance (`1e-10` in the current adapter), so changing
the shared integration tolerances with that family is rejected explicitly.

## Current Limits

The `solve_network` convenience API is one-shot. `PreparedNetworkIvp` borrows
the validated problem and options, reuses the native prepared callbacks, and
can rerun with different initial concentrations over the same time interval.
Constant T or the two linear-heating parameters can be rebound numerically;
arbitrary `Symbolic(T(t))` has a fixed expression and cannot be rebound via
this typed API. Numeric rate coefficients can be rebound without rebuilding
the symbolic graph. Changing the temperature-law shape or interval still
requires a fresh preparation. Sampled output must be finite,
non-negative, time-ordered,
and cover the requested interval. There is no clipping or renormalization.
An opt-in cooperative cancellation token reaches StiffAPI. An aborted run
returns a typed native error and does not replace the last accepted result.
When telemetry is enabled, the latest RST attempt report can be inspected
separately from the accepted solution. It carries a `solve_id`, cancellation
disposition, and optional telemetry, typed callback facts, and at most 32
partial samples. Partial samples are diagnostic, not accepted physical output;
pre-cancelled and failed retries do not inherit the prior attempt's identity or
trajectory. On RST 0.5.2, telemetry and partial-state availability for an
immediate callback failure is family-dependent: Radau and Backward Euler
provide them in the tested case, while BDF and LSODE2 provide only the attempt
identity and typed solve error. Trial time/state and complete callback context
remain unavailable on some native error paths.
Symbolic T(t) is checked at interval endpoints and returned samples; endpoint
checks do not prove positivity at every internal solver trial. Invalid domains
inside the solve are reported by the native backend, not yet enriched with
reaction/time context. Native telemetry is returned only when requested.
When telemetry is enabled, the result also carries a KiThe-owned stage
snapshot for equation assembly, backend preparation, integration, and combined
output postprocessing. Prepared-session stage values are cumulative: setup
stages occur once and integration/postprocessing accumulate over reruns.
Counters contain calls only, Timings contain elapsed durations only, Detailed
contains both, and Off creates no stage snapshot or timing clock reads.

The scalar `solve_1d()` preset is not used for this vector IVP. KiThe checks
model domains at construction, in its explicit RHS probe, and at publication.
Native solver trial states are evaluated by RST; an intermediate native
callback failure may not yet carry stage/model context.
