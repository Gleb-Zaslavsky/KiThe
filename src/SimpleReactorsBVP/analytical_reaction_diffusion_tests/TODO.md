# BVP-01..04 Analytical Reactor Validation

## Scope

Independent correctness coverage for the existing `SimpleReactorTask` equations.
The test-only implementations do not change production physics, transport
closure, or solver algorithms. BVP-02/03 reuse the TW Padé-13 matrix-exponential
action and build operators independently of KiThe's symbolic RHS.

## Checklist

- [x] Freeze the 4-species/5-reaction and 12-species/20-reaction networks,
  component ordering, inlet fractions, and rate constants.
- [x] Independently construct `W` and the `2N x 2N` transport operator; verify
  mass and elemental stoichiometric invariants and estimate `cond2(E_JJ)` by
  singular values without explicitly inverting `E_JJ`.
- [x] Solve the mixed boundary system with pivoted LU and verify its residual.
- [x] Compare production local RHS and symbolic Jacobian against independently
  assembled values, including the interleaved KiThe state ordering.
- [x] Compare Damped/DoublePoints, dense-grid Frozen, and BVP_sci profiles with
  the matrix-exponential oracle for BVP-02a and BVP-02b.
- [x] Measure fixed-mesh convergence for the small network; preserve separate
  tolerances for mass fraction, physical species flux, temperature, and heat flux.
- [x] Compare ExprLegacy and AtomView on the same Damped/DoublePoints BVP-02a
  case and against the same physical oracle.
- [x] Cover unequal species diffusion without composition normalization and
  the zero-kinetics/alpha continuation with the matrix-exponential reference.
- [ ] Release-mode correctness campaign and full BVP regression after the
  debug validation suite is complete.

## BVP-03 Thermochemical Reaction-Diffusion

- [x] Reuse the fixed BVP-02a and BVP-02b reaction fixtures; derive every
  reaction heat from `Q = -S^T h` using the frozen species enthalpies.
- [x] Calibrate `gamma` over `{0, 10, 100, 300, 1e3, 1e4}` on both networks;
  select `gamma=300`, which gives chemical temperature shifts of 28.739 K and
  9.608 K with boundary-system condition estimates `5.838e5` and `1.235e7`.
- [x] Assemble the full block operator, solve mixed thermal/species boundary
  data by pivoted LU, and verify boundary residuals without explicit inversion.
- [x] Compare production RHS and full symbolic Jacobian to the independent
  operator, and finite-difference representative columns.
- [x] Check Damped/DoublePoints, dense Frozen, and adaptive BVP_sci for 4/5 and
  12/20 networks; verify ExprLegacy/AtomView parity on the small network.
- [x] Check internal and dimensional energy first integrals, with term-scaled
  relative defects and max-error coordinates; check heat-release quadrature
  against the species-flux identity and demonstrate trapezoid convergence.
- [x] Check unequal-diffusion energy conservation without normalizing species
  fractions; verify the `gamma=0` production solution reproduces the species
  solution at `gamma=300` and the analytic temperature response is linear in
  `gamma`.
- [ ] Run and record the Release BVP-03 campaign and complete BVP-01/02
  regression before treating these measurements as Release evidence.

Debug evidence is in `../STORY_TESTS.md`. The 12/20 BVP_sci solve takes about
99 seconds in Debug; keep it as one large-network correctness gate.

## Current Measurement Note

For Pe=1 with the prescribed thermal data, Damped at 257 published nodes has
about `7.17e-2 K` temperature and `7.17e-1 W/m2` heat-flux error. Frozen at 1025
nodes reduces these to about `1.80e-2 K` and `1.80e-1 W/m2`; BVP_sci reaches
roughly `1.5e-9 K` and `1.5e-8 W/m2`. These thermal errors match the independent
BVP-01 transport baseline and are reported separately from the much smaller
reactive species errors. The 12-species BVP_sci debug solve is correct but takes
about 105 seconds; keep it as one large-network story gate rather than repeating
that expensive route across every metamorphic case.

## BVP-04 Temperature-Dependent Reaction-Diffusion

- [x] Audit the production Arrhenius contract: `k=A*T^n*exp(-E/(R*T))`,
  `R=8.314`, `E` in J/mol; `T=dT+T_scale*Teta`; first-order rate uses
  `rho*C_i/M_i`; mass source is `M_i*nu_ij*r_j`.
- [x] Freeze SciPy DOP853 fundamental-matrix references for 4/5, 12/20, and
  O3/O2/O/N2 cases. Recompute `A_j` from the BVP-02 `k_ref` values at 330 K;
  check tolerance and max-step convergence, the mixed-boundary solve, its
  condition estimate, and mass/element first integrals.
- [x] Verify Arrhenius values at 298.15, 330, and 361.362 K; compare local
  production RHS and the full symbolic Jacobian, including `dJ'/dTeta`, with
  independent mass-weighted analytic blocks and finite differences.
- [x] Check the isothermal `q_R=0` limit against the constant matrix-exponential
  solution with Damped, Frozen, and BVP_sci.
- [x] Check the small nonautonomous 4/5 case against DOP853 with all three
  solver routes; compare ExprLegacy and AtomView profiles.
- [x] Check the unequal-molar-mass O3/O2/O/N2 case against DOP853 with all
  three solver routes; assert molecular masses, mass balance, O/N atom balance,
  and solution-level first integrals.
- [ ] Complete the 12/20 BVP_sci production-profile gate. A Debug attempt with
  `max_nodes=128` failed with `MaxNodesExceeded`; retrying with 512 nodes did
  not complete within four minutes and was stopped. The test remains ignored
  for ordinary Debug runs and must be run in Release or recalibrated without
  weakening the frozen-reference acceptance.
- [ ] Run the Release BVP-04 campaign plus BVP-01/02/03 regression.

The dimensional temperature reference uses the physical outlet heat flux:
`T(z)=T0 + q_phys,R*L/lambda*(exp(Pe_q*(z-1))-exp(-Pe_q))/Pe_q`.
In the dimensionless equation the boundary value is instead
`q_hat_R=q_phys,R*L/T_scale`; do not insert the physical `1000 W/m2` directly
as the internal `q` value.
