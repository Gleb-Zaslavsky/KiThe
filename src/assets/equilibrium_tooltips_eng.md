## tab.setup
Define the thermodynamic problem, inventory, phases, and initial amounts.

## tab.phase_control
Configure phase activation, deactivation, TPD, hysteresis, and lifecycle budgets.

## tab.libraries
Choose lookup sources and inspect offline and fallback provenance policy.

## action.run
Resolve the request and run the selected equilibrium workflow.

## action.cancel
Cancel the running calculation without publishing a partial result.

## mode.pt
Solve equilibrium at fixed pressure and temperature. The selected mode is
highlighted and determines which canonical facade request is prepared.

## mode.ph
Solve equilibrium at fixed pressure and total enthalpy. The selected mode is
highlighted and determines which canonical facade request is prepared.

## temperature.point
Solve one fixed temperature in kelvin.

## temperature.range
Evaluate an ordered temperature grid with continuation between accepted points.

## output.pchip
Resample a P,T range with shape-preserving PCHIP for display only.

## diagnostics.timing
Collect stage and per-point timing evidence for the calculation.

## action.validate
Validate the editable document without starting a solve.

## action.prepare
Resolve data and build the canonical equilibrium request.

## action.embedded_plot
Open a plot from the accepted snapshot without solving again.

## action.kithe_plot
Open the KiThePlot view from the accepted snapshot.

## field.target_enthalpy
Target total enthalpy in joules for the P,H constraint.

## field.moles
Initial amount in moles; it must be finite and non-negative.

## field.element
Chemical element symbol used by element-based candidate search.

## field.phase_epsilon
Small phase amount threshold used by lifecycle decisions.

## field.dg_create
TPD threshold below which an absent phase may be created.

## field.dg_keep
TPD threshold controlling retention of an active phase.

## field.phase_iterations
Maximum number of outer phase-lifecycle iterations.

## field.trace_floor
Absolute positive trace-species seed in moles.

## field.trace_fraction
Trace seed as a fraction of the largest initial amount.

## field.trace_minimum_floor
Lower bound applied to a relative trace-species seed.

## field.tolerance
Optional nonlinear acceptance tolerance override.

## field.max_iterations
Optional per-backend nonlinear iteration limit.

## field.cascade_attempts
Maximum number of fallback backend attempts.

## field.cascade_iterations
Iteration budget for each backend attempt.

## field.cascade_total
Total iteration budget across the fallback cascade.

## field.display_points
Number of display-only interpolation points for P,T ranges.

## inventory.explicit
Enter phase-qualified substances and their physical initial amounts directly.

## inventory.elements
Use the requested elemental inventory to search permitted libraries for
candidates. Candidate search does not create synthetic thermodynamic species.

## inventory.state
Select gas, liquid, or solid semantics for this phase.

## inventory.model
Choose the thermodynamic model used by this phase.

## lookup.priority
Libraries are tried in this ordered priority list.

## lookup.priority_entry
A library name in the ordered lookup priority list.

## lookup.permitted
Restrict candidate lookup to this closed set of libraries.

## lookup.permitted_entry
A library name allowed by the closed candidate lookup set.

## lookup.nist
Allow an explicit NIST online fallback when local lookup is insufficient.

## lookup.default
Use the repository's production lookup policy.

## lookup.explicit
Configure the ordered and permitted libraries explicitly.

## lookup.load_catalog
Load local library names without resolving substances or starting a solve.

## lookup.catalog
Choose a locally discovered library name for explicit lookup policies.

## candidate.states
Restrict candidate search to the selected physical states.

## candidate.max
Limit the number of candidates returned by the preview.

## candidate.state.gas
Include gas-phase records in the candidate search.

## candidate.state.liquid
Include liquid-phase records in the candidate search.

## candidate.state.solid
Include solid-phase records in the candidate search.

## candidate.state.condensed
Include condensed-state records in the candidate search.

## phase.initial_positive
Start lifecycle control from phases with positive initial inventory.

## phase.initial_all
Also test declared zero-inventory phases as appearance candidates.

## phase.initial_set
Choose whether the initial active phase set follows inventory or all declared candidates.

## solver.production
Use the production nonlinear backend cascade and fallback policy.

## solver.backend
Use one concrete nonlinear backend for diagnosis or a custom cascade.

## solver.custom
Create an editable ordered nonlinear fallback cascade.

## solver.trace_absolute
Seed every trace species with this absolute positive mole floor.

## solver.trace_seed
Choose how trace-species initial seeds are scaled.

## solver.trace_relative
Scale trace seeds relative to the largest initial mole amount.

## diagnostics.phase_off
Do not retain phase-lifecycle events.

## diagnostics.phase_summary
Retain a compact summary of phase-lifecycle decisions.

## diagnostics.phase_lifecycle
Retain activation, deactivation, and accepted lifecycle events.

## diagnostics.phase_detailed
Retain detailed TPD and stability evidence; this uses more memory.

## diagnostics.phase_trace
Select how much phase-lifecycle evidence is retained.

## diagnostics.range_trace
Select which temperature-range lifecycle points are retained.

## preset.ideal_gas
Replace the inventory with a small valid ideal-gas example.

## diagnostics.range_endpoints
Retain lifecycle evidence for the first and last range points.

## diagnostics.range_transitions
Retain only range points with phase transitions.

## diagnostics.range_every_point
Retain lifecycle evidence for every range point; this uses more memory.

## tab.numerics
Choose the nonlinear backend cascade and numerical safeguards.

## tab.output
Configure retained diagnostics, result presentation, and plotting.

## tab.results
Inspect the immutable accepted result snapshot and its evidence.

## help.language
Select the language used by equilibrium help and hover explanations.

## language.english
Show calculator help and hover explanations in English.

## language.russian
Show calculator help and hover explanations in Russian.

## candidate.policy
Control element matching and the size of the candidate preview.

## candidate.matching
Choose subset matching or require exactly the requested element set.

## candidate.bounds
Limit candidate records to this thermochemistry temperature interval.

## candidate.preview
Resolve and display candidates without starting equilibrium solving.

## candidate.cancel
Cancel candidate lookup; no partial preview is published.

## candidate.target
Choose the phase receiving selected candidate assignments.

## list.up
Move this row earlier in the ordered list.

## list.down
Move this row later in the ordered list.

## list.remove
Remove this row from the editable configuration.

## list.add_phase
Add another phase candidate to the problem.

## list.add_component
Add a component to this phase.

## list.assign
Assign the selected catalog candidate to the target phase.

## list.unassign
Remove the candidate assignment without changing catalog data.

## list.add_backend
Append a backend to the ordered fallback cascade.

## list.add_library
Append the selected library to this policy list.

## list.add_element
Add an element symbol to the candidate-search query.

## phase.fixed
Keep the declared phase set fixed during the solve.

## phase.bounded
Allow activation and deactivation using TPD and hysteresis.

## solver.ph_route
Choose Auto, monolithic, or nested temperature handling for P,H.

## solver.ph_auto
Let the production policy choose the P,H route and fallbacks.

## solver.ph_monolithic
Solve composition and temperature together in one nonlinear system.

## solver.ph_nested
Solve P,T inner states while an outer search determines temperature.

## solver.budget
Override the default fallback cascade iteration budget.

## solver.scaling
Enable numerical scaling for better conditioning of the solve.

## solver.trace_override
Override the default trace-species seed policy.

## output.basis
Choose the physical quantity displayed and plotted.

## output.plot_target
Choose where an accepted result is plotted.

## output.plot_none
Do not open a plot after an accepted solve.

## output.plot_embedded
Use the calculator's embedded plot view.

## output.plot_kithe_plot
Send accepted data to the KiThePlot editor.

## output.plot_both
Open both plot targets from the same accepted snapshot.

## output.scale
Choose the vertical display scale for plots.

## output.table_density
Choose how much derived information each visible result row shows.

## output.table_detailed
Show moles and mole fractions in the result table.

## output.table_compact
Show component moles only in the result table.

## output.visible_series
Show or hide this series in presentation and plots only.

## output.interpolation
Choose linear or logarithmic interpolation for display only.

## output.clamp
Keep interpolated display values inside the solved range.

## diagnostics.keq
Use the independent equilibrium-constant validator when applicable.

## diagnostics.keq_off
Disable independent equilibrium-constant validation.

## diagnostics.keq_when_applicable
Run K_eq validation only for compatible small systems.

## diagnostics.keq_required
Reject a result when applicable K_eq validation cannot be completed.

## state.gas
Gas phase; ideal-gas composition and pressure activities apply.

## state.liquid
Liquid phase candidate.

## state.solid
Solid phase candidate.

## model.ideal_gas
Ideal gas thermodynamic model.

## model.ideal_solution
Ideal solution model for a mixed phase.

## model.pure_condensed
Pure condensed-phase model.

## matching.subset
Allow records containing the requested elements plus additional elements.

## matching.exact
Require the record element set to equal the requested set.

## basis.moles
Plot component amounts in moles.

## basis.fractions
Plot component mole fractions.

## basis.phase_totals
Plot total amount per phase.

## scale.linear
Display values on their ordinary linear scale.

## scale.log
Display positive values on a base-10 logarithmic scale.

## interpolation.linear
Interpolate display values linearly.

## interpolation.log
Interpolate positive display values in log space.

## field.pressure
System pressure in pascals. It must be finite and positive; changing it invalidates preparation.

## field.reference_pressure
Gas activity reference pressure in pascals. It is distinct from system pressure and must match the data convention.

## field.phase_id
Stable phase identifier used to distinguish this phase in results and plots.

## field.substance
Thermochemical component identifier resolved under the selected library policy.

## field.display_cutoff
Presentation cutoff for small mole fractions. It filters visible rows without changing solved data.
