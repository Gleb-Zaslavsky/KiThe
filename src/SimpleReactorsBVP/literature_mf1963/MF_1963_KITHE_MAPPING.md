# MF-1963: KiThe Input Mapping and Current Integration Evidence

This is the implementation-facing companion to `MF_1963_REFERENCE.md`. It
records how the frozen dimensional model is expressed by the existing
`SimpleReactorTask` API. It does not change production equations or solver
defaults.

## Model encoding

The two formal species are named `A` and `B`, with identical elemental
composition `C2H4O2` and molar mass `0.060052 kg/mol`. This is an isomolar
one-to-one conversion, so KiThe mass fractions equal the reference conversion
and unconverted reactant fraction: `Y_B=eta`, `Y_A=1-eta`. The reaction is
`A=>B`. In `FastElemReact`, `n` is the temperature exponent; the kinetic order
comes from the reactant stoichiometry, so `A=>B` is first order with `n=0`.

`FastElemReact.Q` and the parallel `thermal_effects` entry are molar heats in
J/mol. The conversion used is
`Q_mol = 1.5e6 J/kg * 0.060052 kg/mol = 90078 J/mol`. The reaction's
mass-specific heat release is therefore recovered by dividing by the common
molar mass. The inlet has `A=1`, `B=0`, and `T=300 K`; the outlet heat and
species fluxes are zero.

The affine temperature transform is `T=300 + 15*Teta`, so
`Teta=theta`. Each case gets its own physical length from MF-00 metadata:

```text
L = lambda*xi_max/(m*cp)
Pe_q = m*cp*L/lambda = xi_max
xi = Pe_q*z,  z=x/L
```

The `DensityModel::IdealGasFixedComposition` default is retained, with
`Tm=T0=300 K` and pressure `P=rho0*R_G*Tm/M`, making the reference density
`rho0=1 kg/m^3`. The production transport code computes
`rhoD = D_input*rho_298*sqrt(Tm/298.15)`, where
`rho_298=M*P/(R_G*298.15)`. Therefore the API input is inverted as

```text
D_input = (lambda/cp)/(rho_298*sqrt(Tm/298.15))
```

The preparation test confirms both species recover `rhoD=1e-4 kg/(m s)` and
`Le=1`.

The profile seed maps the reduced state into the actual six KiThe unknowns
`[Teta, q, C0, J0, C1, J1]`. `q` is seeded as
`lambda*dTeta/dz = lambda*Pe_q*p`; species flux states are
`J_A=(rhoD)*dY_A/dz` and `J_B=(rhoD)*dY_B/dz`. Postprocessing converts them
back with `q_physical=q*T_scale/L` and `J_physical=J/L`. MF-00 seeds are not
the only starting point: MF-HIGH also has an independently constructed smooth
profile seed.

## Gas-constant convention

MF-00 uses `R=8.31446261815324 J/(mol K)`, while KiThe's BVP model uses
`R_G=8.314`. MF-01 chooses `E=R_G*T0/beta` and reconstructs `k0` with the same
`R_G`, so KiThe recovers the requested `beta=0.05`, `gamma=0.01`, and each
requested `delta`. Evaluating these same inputs under MF-00's R gives
`beta=0.0500027822`, `gamma=0.0100005564`, and
`delta_MF00/delta_KiThe=1.00105772`. This approximately `0.106%` delta
difference is recorded rather than hidden; the comparison is against the
actual MF-00 fixture, not a claim of exact constant parity.

## Current routes and findings

The detailed current route/configuration/outcome matrix is **MF-07** in
[`STORY_TESTS.md`](../STORY_TESTS.md#mf-07-highlow-production-capability-and-rst-parity-summary-2026-11-12).
It supersedes earlier initial MF-01/MF-02 solver-selection observations below
this heading. In particular, the old claims that no MF-LOW production route
passes, that MF-LOW BVP_sci fails, and that the tested Damped strategies all
select the cold branch are no longer current: with the explicit clustered mesh
and revised settings, MF-LOW BVP_sci and Damped/DoublePoints pass the physical
gate. Solver convergence and physical-reference acceptance remain separate
columns in the current matrix.

## RST 0.4.35 Explicit-Mesh Update

The older fixed-mesh limitation described in MF-04 is superseded. KiThe now
offers `BVPSolver::solve_with_bvp_api_on_mesh`, which accepts a solver-coordinate
mesh and a full state-major `state_count x node_count` profile including both
boundary nodes. The method validates mesh ordering/domain, profile shape and
finiteness, and consistency of boundary-aligned seed values before the normal
RST handoff. It keeps ordinary transactional publication and dimensional
postprocessing. Damped starts from this mesh and may refine; CorrectionNewton
keeps it fixed; BVP_sci can adapt its collocation mesh.

The current production MF-LOW snapshot is in `STORY_TESTS.md` under MF-05 and
the consolidated HIGH/LOW conclusion is MF-07.
BVP_sci and Damped/DoublePoints with four refinements pass the physical
reference-profile gate. CorrectionNewton/BandedAuto, Frozen with
`at_high_morm=1e-6`, and Damped TwoPoint/Pearson/Grcar-Smooke return converged
discrete roots but miss that gate. Frozen's LOW profile matches the measured
CorrectionNewton alternate branch; the remaining three Damped policies produce
less distant but still rejected profiles. This distinction is asserted through
the production facade, not inferred from solver status alone.

For MF-HIGH, the 301-node uniform production matrix passes for
CorrectionNewton/BandedAuto, Damped TwoPoint/Pearson/Grcar-Smooke/DoublePoints,
and Frozen default/`at_high_morm=1e-6`/`every_m=1`. The Frozen policies agree
with the fixed-mesh CorrectionNewton profile within `1.1e-8`. These findings
are specific to the model, seed, mesh, and settings. The high/low solver-option
transfer is now substantially covered; production cross-seed coverage,
release-campaign confirmation, and the remaining duplicate audit are still
open.

## Source instruction status

The supplied MF-01 pasted document ends mid-equation in section 8, before its
remaining acceptance criteria and route matrix. The tests and mappings above
are based on the complete available requirements and the current code. Do not
infer or mark the missing portions complete; reconcile them if the complete
instruction becomes available.
