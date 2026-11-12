# Merzhanov-Filonenko 1963: Independent BVP Reference

Status: MF-00 mathematical reference. This document and its fixtures are
generated independently with NumPy/SciPy. No KiThe or RustedSciThe solver is
used here; this is not yet a KiThe validation result.

## Source and model scope

Reference: A. G. Merzhanov and A. K. Filonenko, “О тепловом самовоспламенении
гомогенной газовой смеси в потоке”, Доклады АН СССР, 1963, 152(1), 143–146.

The selected model is a steady, one-dimensional, constant-pressure flow with
constant mass flux and a single irreversible first-order reaction. The original
Arrhenius law is retained; no Frank-Kamenetskii approximation is made. Specific
heat `cp`, thermal conductivity `lambda`, and `rho*D` are constant, with
`Le = 1` and `rho*D = lambda/cp`. Density changes with temperature as
`rho(T) = rho0*T0/T`; composition effects on density are omitted.

Dimensional equations, using conversion `eta` (zero at the inlet and one after
complete conversion), are

```text
lambda * d2T/dx2 - m * cp * dT/dx + Q * W = 0
rho*D * d2eta/dx2 - m * deta/dx + W = 0
W = rho(T) * k0 * exp(-E/(R*T)) * (1 - eta)
```

The reaction source is positive in the conversion equation. If reactant mass
fraction is used instead, `Y_A = 1 - eta` and its source has the opposite sign.
Boundary conditions on the semi-infinite domain are `T(0)=T0`, `eta(0)=0`,
`dT/dx -> 0`, and `deta/dx -> 0` downstream. The numerical reference uses a
large finite interval and imposes the two inlet values and zero downstream
derivatives.

Because `Le=1` and the transport coefficients are constant, the equations
admit the independent first-integral check
`eta = cp*(T-T0)/Q`.

## Dimensionless system

Define

```text
theta = E/(R*T0^2) * (T - T0)
xi    = (m*cp/lambda) * x
beta  = R*T0/E
gamma = cp*R*T0^2/(Q*E)
delta = Q*lambda*rho0*E/(cp^2*m^2*R*T0^2) * k(T0)
k(T0) = k0*exp(-E/(R*T0))
```

Then `T/T0 = 1 + beta*theta`, `eta = gamma*theta`, and the first-order system
is

```text
theta' = p
p'     = p - S(theta)
S(theta) = delta * exp(theta/(1 + beta*theta))
                   * (1 - gamma*theta)/(1 + beta*theta)
```

The finite-domain boundary conditions are `theta(0)=0` and `p(xi_max)=0`.
The outlet temperature is an output, not a prescribed Dirichlet boundary
condition. In particular, the solver does not impose `theta(xi_max)=1/gamma`.

The source derivative used by the analytic collocation Jacobian is obtained
from `g=theta/(1+beta*theta)`:

```text
S' = delta*exp(g) * (
       (1-gamma*theta)/(1+beta*theta)^3
       - gamma/(1+beta*theta)
       - beta*(1-gamma*theta)/(1+beta*theta)^2
     )
```

The implementation evaluates this expression directly, including near full
conversion; it does not divide by `1-gamma*theta`.

## Parameters and dimensional mapping

The reference uses `beta=0.05`, `gamma=0.01`, reaction order `n=1`, and the
following dimensional baseline:

| Parameter | Value |
| --- | ---: |
| `T0` | 300 K |
| `cp` | 1000 J/(kg K) |
| `rho0` | 1 kg/m^3 |
| `lambda` | 0.1 W/(m K) |
| mass flux `m` | 0.01 kg/(m^2 s) |
| `E` | `R*T0/beta` = about 49.887 kJ/mol |
| `Q` | `cp*T0*beta/gamma` = 1.5 MJ/kg |
| `rho*D` | `lambda/cp` = 1e-4 kg/(m s) |
| adiabatic `T` | `T0 + Q/cp` = 1800 K |

The dimensional pre-exponential is reconstructed per case from
`k(T0) = delta*cp*gamma*m^2/(lambda*rho0)` and
`k0 = k(T0)*exp(E/(R*T0))`. The exact values and floating-point environment are
written to `fixtures/metadata.json` by the generator.

The dimensional inlet heat flux is reported using
`q_in = lambda*dT/dx|0 = T_scale*m*cp*p(0)`. The physical axial coordinate is
`x = lambda/(m*cp)*xi`.

For a later isomolar KiThe reaction mapping, `FastElemReact.Q` is molar heat in
J/mol. Use `Q_mol = Q_mass*M`, where `M` is the common molecular mass in kg/mol
of the isomolar reactant/product pair. The model reference itself remains on a
mass basis. The suggested future axial mapping is
`x = lambda/(m*cp)*xi`, `T = T0 + (R*T0^2/E)*theta`, and `Y_A=1-gamma*theta`.
The density mapping is `rho(T)=rho0*T0/T`. A future implementation must verify
KiThe's exact sign and nondimensionalization conventions rather than assume
these formulas are wired identically.

## Frozen scenarios

The two cases are deliberately separated by their resolved profiles:

| Case | `delta` | Intended regime |
| --- | ---: | --- |
| `MF-LOW` | `1e-4` | Long induction and weak inlet heat flux; the reaction front lies far downstream. |
| `MF-HIGH` | `0.1` | Intense, compact source/front near the inlet and much larger inlet heat flux. |

Selection is supported by computed front location, inlet slope, source maximum,
and integrated source, not only by the parameter values. These are two measured
solutions, not a critical-parameter study and not a proof of uniqueness.

## Numerical method and acceptance checks

`generate_reference.py` uses SciPy `solve_bvp` with the analytic Jacobian of
the first-order RHS. A stable-manifold trajectory near the adiabatic equilibrium
is used only as a mesh/initial-profile seed so that the high-activation
collocation problem starts on the physical branch. SciPy's collocation solve
then enforces the finite-domain BVP. Both the adaptive mesh and the physical
domain are varied independently; a nonzero solver status is a hard failure.

The local checks require finite values, `1+beta*theta > 0`,
`0 <= gamma*theta <= 1` up to numerical tolerance, inlet/outlet boundary
conditions, and agreement between the independently quadrature-computed source
integral and the ODE identity

```text
integral(S dxi) = theta(xi_max) - theta(0) - p(xi_max) + p(0).
```

The integral is evaluated with adaptive quadrature over SciPy's continuous
collocation solution; it is not substituted by the identity under test. The
low case also reports a first-order convective approximation as a qualitative
cross-check. Literature figures are qualitative context only; no digitized
curve or machine-precision match is claimed.

## Reproduction

From the repository root:

```powershell
uv run --with numpy --with scipy --with matplotlib python src/SimpleReactorsBVP/literature_mf1963/generate_reference.py
uv run --with numpy --with scipy --with matplotlib python -m unittest discover -s src/SimpleReactorsBVP/literature_mf1963 -p 'test_*.py' -v
```

The generator writes deterministic 17-digit profile CSVs, profile SVGs, and
solver metadata under `fixtures/`. The low-case convective approximation
compares the first-order prediction for the location where `theta=90` against
the full BVP; it is a qualitative cross-check, not a replacement solution.
The tests validate the equations and frozen artifacts;
they do not invoke a KiThe solve. Numerical values, SciPy version, mesh/domain
studies, and final case diagnostics are appended here after generation.
