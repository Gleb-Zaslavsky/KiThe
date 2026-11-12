# Gas-Phase Reactor BVP Model

This document records the physical balances, solver coordinates, boundary
conversion, and publication rules implemented by the gas-phase BVP route. The
equation assembly lives in `createBVP.rs`; scale conversion and result
publication live in `SimpleReactorBVP.rs` and `SimpleReactorBVP3.rs`.

## Coordinates And State

The physical axial coordinate is `x` in metres and `L` is the characteristic
reactor length. The BVP solver uses

```text
z = x / L,        x = L z.
```

The task's `t0` and `t_end` are therefore bounds in `z`, despite the legacy
field names. The solver argument label is caller-selected; the numerical mesh
is still dimensionless until postprocessing.

The thermal state is

```text
Teta = (T - dT) / T_scale,        T = dT + T_scale Teta,
```

where `dT` is a temperature offset and `T_scale` is a positive temperature
scale, both in kelvin. They are distinct parameters and need not be equal.

Each species state `C_i` is a mass fraction, not a molar concentration:

```text
C_i = w_i,        0 <= C_i <= 1,        sum_i C_i = 1 at the inlet.
```

For reaction kinetics, the code converts this to molar concentration using
`c_i = rho(T) C_i / M_i`, where `rho(T)` follows the selected density closure
and `M_i` is the molar mass in `kg/mol`.

## Density And Transport Closures

`DensityModel::IdealGasFixedComposition` is the default. It assumes constant
pressure and a fixed reference mixture molar mass:

```text
rho_ref = P M_ref / (R Tm),
rho(T) = rho_ref Tm / T,
T = dT + T_scale Teta.
```

`M_ref` is the task's reference mean molar mass; it is not recomputed from the
local composition. The alternative `DensityModel::Constant` explicitly freezes
density at `rho_ref`.

Both density modes retain the existing transport closure `K_i = rho D_i =
constant`, with unchanged `D_ro_map` and `Pe_D`. Consequently, in the
temperature-dependent mode `D_i(T) = K_i / rho(T)` is implied by the closure;
KiThe does not separately evaluate a temperature-dependent diffusion model.
Thermal conductivity, heat capacity, reaction heats, `Pe_q`, and the boundary
condition contract are unchanged.

The solver's flux states are coordinate-scaled states, not physical fluxes:

```text
q_hat = Lambda d(Teta)/dz,
J_hat_i = rho D_i d(C_i)/dz.
```

With the implementation's positive-gradient convention, their dimensional
counterparts are

```text
q_physical = (T_scale / L) q_hat       [W/m^2],
J_physical_i = J_hat_i / L             [kg/(m^2 s)].
```

These definitions do not include the conventional leading minus sign sometimes
used for Fourier or Fick fluxes. A positive configured flux follows the
gradient convention above.

## Reaction Sources

The stoichiometric coefficient `nu_ij` is signed (reactants negative, products
positive). For reaction rate `r_j` in `mol/(m^3 s)`, the mass source for species
`i` is

```text
S_i = M_i sum_j(nu_ij r_j)              [kg/(m^3 s)],
r_j = k_j(T) product_i[(rho C_i/M_i)^alpha_ij].
```

Rate-constant units depend on the reaction order and must make `r_j` have the
stated volumetric molar-rate units. Reaction heat `Q_j` is supplied in
`J/mol`; the assembled volumetric heat source is `Q_vol = sum_j(Q_j r_j)`.

## Dimensional Balances

With the configured constant transport products `K_i = rho D_i`, constant
`Lambda`, `Cp`, and mass flux `m`, the dimensional balances represented by the
assembled RHS are

```text
d/dx(K_i dC_i/dx) - m dC_i/dx + S_i = 0,
d/dx(Lambda dT/dx) - m Cp dT/dx + Q_vol = 0.
```

For the default variable-density model, `rho(T)` affects the reaction
concentrations and sources, while `K_i` remains constant. The equations assume
steady one-dimensional transport and no additional source or loss terms beyond
the configured reactions.

## Solver Equations

Using `z=x/L`, `Teta=(T-dT)/T_scale`, and the flux-state definitions above,
the species equations assembled by KiThe are

```text
dC_i/dz = J_hat_i / K_i,
dJ_hat_i/dz = Pe_D_i J_hat_i - L^2 S_i,
Pe_D_i = m L / K_i.
```

The thermal equations assembled by KiThe are

```text
dTeta/dz = q_hat / Lambda,
dq_hat/dz = Pe_q q_hat - (L^2/T_scale) Q_vol,
Pe_q = m Cp L / Lambda.
```

The source signs above describe the current production equations and the
configured signed stoichiometry/heat convention. They are a faithful statement
of the implementation, not an independent validation of every physical sign
convention.

## Boundary Conditions

The task-facing boundary conditions are physical inputs where a dimensional
quantity exists:

| Task field | Meaning and units | Solver value |
| --- | --- | --- |
| Species name, e.g. `HMX` | Inlet mass fraction, dimensionless | `C_i = w_i` |
| `T` | Inlet temperature, K | `Teta = (T - dT) / T_scale` |
| `q` | Outlet conductive heat flux, W/m^2, using the positive-gradient convention above | `q_hat = q L / T_scale` |
| `<species>_flux` | Outlet species mass flux, kg/(m^2 s), using the positive-gradient convention above | `J_hat_i = J_i L` |

Absent `q` and `<species>_flux` fields retain the existing small solver-space
near-zero defaults (`1e-10`) for the adiabatic/closed endpoint. Explicit values
are checked after conversion as well as before it, so overflow is rejected.

Solver bounds and initial guesses (`C`, `J`, `Teta`, `q`) are numerical
solver-space values; they are not task boundary inputs and are not converted.

## Publication And Units

The standard task solve applies dimensional publication by default
(`postprocessing.return_to_dimension: true`):

```text
z       -> x = L z,
Teta    -> T = dT + T_scale Teta,
q_hat   -> q_physical = q_hat T_scale / L,
J_hat_i -> J_physical_i = J_hat_i / L,
C_i     -> unchanged mass fraction.
```

With `return_to_dimension: false`, the stored and displayed profile remains in
solver coordinates: dimensionless `z` and `Teta`, mass fractions `C_i`, and
coordinate-scaled flux states `q_hat` and `J_hat_i`. This option is explicit;
it is not the default physical-output route. Molar fractions and molar
concentrations are separate derived quantities and are not substituted for the
published mass-fraction state automatically.

The GUI uses the solution snapshot produced by the selected postprocessing
setting and must not apply the dimensional transform a second time. Its plot
labels distinguish `Temperature (K)`, `Teta (dimensionless)`, physical `q/J`
units, solver-coordinate `q_hat/J_hat`, and species mass fractions.

## Scale Validation

`dT`, `L`, and `T_scale` must each be finite and strictly positive. `dT` is the
temperature offset; `T_scale` is the divisor/multiplier in the temperature
transform and in the thermal source coefficient. A zero, negative, NaN, or
infinite value does not define a valid scaled problem.
