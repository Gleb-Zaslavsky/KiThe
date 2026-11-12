# Condensed Reactor IVP Model

This note records the physical assumptions and the exact dimensionless system
assembled by the condensed-phase IVP route. It is the reference for the RHS
builder, fixtures, and contract tests.

## Coordinates And State

The physical coordinate is `x` in metres and the characteristic length is `L`.
The IVP solver integrates the dimensionless coordinate

```text
z = x / L,        x = L z.
```

Accordingly, integration bounds currently stored under the compatibility names
`x0` and `x_bound` are values of `z`. The raw solver argument is named `z`;
published profiles convert it to physical `x` by multiplying by `L`.

The species states `Y_i` are mass fractions. The thermal state is
`Teta = (T - dT) / T_scale`. The solver heat-flux state is
`qhat = q_physical L / T_scale`, so the user-facing initial heat flux in
`W/m^2` is converted to `qhat` before integration and converted back when
publishing results.

Physical inputs are constant density `rho` in `kg/m^3`, positive mass flux `m`
in `kg/(m^2 s)` directed toward increasing `x`, heat capacity `Cp` in
`J/(kg K)`, conductivity `Lambda` in `W/(m K)`, and length `L` in metres.

## Reaction Sources

The stoichiometric analyzer uses signed coefficients `nu_ij`: reactants are
negative and products positive. A forward rate `r_j` is non-negative and has
units `mol/(m^3 s)`. For mass fraction `Y_i` and molar mass `M_i` in `kg/mol`,
the species mass source is

```text
omega_i = M_i sum_j(nu_ij r_j)       [kg/(m^3 s)].
```

For a steady one-dimensional convective balance with no species diffusion,

```text
m dY_i/dx = omega_i.
```

Applying `z = x/L` gives the dimensionless-coordinate equation used by the IVP:

```text
dY_i/dz = (L/m) omega_i
        = (L/m) M_i sum_j(nu_ij r_j).
```

The factor `L/m` converts the volumetric source to a change in mass fraction
over the dimensionless reactor length. The sign follows directly from the
signed stoichiometry: for `A -> B`, `dY_A/dz < 0` and `dY_B/dz > 0`. For a
mass-balanced reaction, `sum_i M_i nu_ij = 0`, so the reaction source preserves
the sum of species mass fractions. Elemental conservation follows when the
reaction stoichiometry conserves each element.

The kinetic expression currently assembled for each reaction is

```text
r_j = k_j(T) product_i[(rho Y_i / M_i) ^ alpha_ij],
```

where `alpha_ij` comes from the reactant-order matrix. The rate-constant units
must match the total reaction order so that `r_j` has the stated volumetric
molar-rate units.

For the first-order TW validation mechanism, all species have the same molar
mass and each reaction consumes one species to produce one species. Therefore
the independently assembled linear mass-fraction system is

```text
dY_i/dz = (L * rho / m) * sum_j(nu_ij * k_j(T) * Y_reactant(j)).
```

Here `nu_ij` is the signed product-minus-reactant coefficient. This equation is
the specialization of the production balance above; it is not the older
concentration-coordinate expression `-rho * L^2 * D_M * S * ...`. The latter
does not represent the current `Y_i` state or its `z=x/L` scaling.

## Arrhenius Kinetics And TW-03

The elementary kinetic backend evaluates

```text
k_j(T) = A_j * T^n_j * exp(-E_j / (8.314 * T)).
```

For TW-03, `n_j=0`, `T` is the absolute kelvin value substituted from
`T=300+100*Teta`, `E_j` is supplied in `J/mol`, and the actual KiThe
`ElementaryStruct` uses `R=8.314 J/(mol K)`. There is no hidden conversion from
kJ/mol. With a reference speed `k_j,ref` at 400 K, the fixture sets

```text
A_j = k_j,ref * exp(E_j / (8.314 * 400)),
k_j(T) = k_j,ref * exp((E_j/8.314) * (1/400 - 1/T)).
```

Thus changing `E_j` leaves the speed at 400 K unchanged. TW-03 uses the same
12 pseudo-species and 20 first-order reactions as TW-01, zero reaction heats,
and an analytical temperature field independent of chemistry. The supplied
draft's chemical generator `-rho*L^2*...` was not copied: KiThe's current
mass-fraction RHS contains `L/m`, and with equal molar masses reduces to the
expression above.

The current TW-01 physical parameters are retained (`rho=1000`, `m=1`,
`Lambda=10000`, `L=0.01`, `Cp=1000`). They give `Pe_q=0.001`. For the
variable-temperature TW-03 cases the internal initial state is `qhat=5000`,
which gives `qhat/Lambda=0.5` and reproduces the requested analytic temperature
profile from 400 K to about 450.025 K. The task builder accepts dimensional
initial `q_physical`; the test supplies its equivalent
`q_physical=qhat*T_scale/L` rather than passing internal `qhat` as if it were
dimensional. This preserves the TW-01 chemical scale; using the attachment's
alternative `m=1e-4`, `Lambda=1` would change the species residence scale by
four orders of magnitude.

TW-03a uses a constant generator and matrix exponential. TW-03b has equal
activation energies, so `K(z)=g(z)K_ref`; an independently integrated scalar
`sigma=integral(g dz)` reduces it to `exp(sigma*K_ref)Y0`. TW-03c uses different
activation energies, checks a nonzero commutator, and uses a separate
fourth-order Magnus reference with a measured refinement study. Frozen profile
data and provenance are in `combustion_validation/fixtures/tw_03_*`.

TW-04 activates the existing thermal coupling and uses the same 12 species, 20 reactions,
enthalpies, activation energies, reference rates, and physical parameters. It sets
`Q_j(gamma) = gamma * (h_from - h_to)`, which is `-gamma * S^T h` for the production signed
stoichiometry. The calibrated fixtures use `gamma = 0`, `1e3`, `1e5`, and `5e5`; the latter is
accepted only because the independent solver comparison, positive mass fractions, and conservation
checks pass. The measured maximum temperature departures from the zero-feedback profile are
approximately `0`, `0.0187`, `2.06`, and `21.9 K`; the strong case departs by about `18.5 K` from
the small-gamma linear extrapolation. Full RHS and trajectory references are in
`combustion_validation/fixtures/tw_04_*`, generated by the independent NumPy/SciPy implementation
under `combustion_validation/external_reference/tw_04`. This validates the declared mathematical
equations, but does not settle the separate physical sign audit for heat flow.

## Thermal Pair

The initial heat flux is supplied as physical `q_physical` in `W/m^2`; the
solver state is `qhat = q_physical L/T_scale`. The existing scaled thermal pair
is

```text
dTeta/dz = qhat / Lambda.
```

For reaction heat `Q_j` in `J/mol`, define volumetric heat release

```text
Qvol = sum_j(Q_j r_j)                 [W/m^3],
Pe_q = m Cp L / Lambda.
```

The currently implemented thermal balance is

```text
dqhat/dz = Pe_q qhat - (L^2/T_scale) Qvol.
```

The factors are dimensionally consistent with the scaled heat-flux coordinate
and volumetric heat source. This records the current RHS, not an independent
physical derivation of its sign. In particular, the relation between positive
`q_physical`, positive coordinate direction, and the configured exothermic
convention `Q_j > 0` still needs its own physical sign audit before we claim a
heat-flow direction.

## Invariants And Test Scope

- A closed, mass-balanced reaction network preserves total mass fraction.
- A closed, element-balanced reaction network preserves each elemental mass
  inventory formed from the species mass fractions.
- The minimal `A -> B` test checks both source signs and the exact `L/m`
  coefficient. The TW analytical fixtures independently implement the same
  dimensionless species generator and thermal pair.
- These equations assume steady one-dimensional convection, constant density
  and mass flux, and no species diffusion. Changing those physical assumptions
  requires deriving a different balance rather than adjusting numerical
  tolerances.
