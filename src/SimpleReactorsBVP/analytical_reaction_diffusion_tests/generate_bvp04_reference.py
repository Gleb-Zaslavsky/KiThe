"""Generate independent nonautonomous BVP-04 checkpoints with SciPy DOP853.

Run from the repository root:
    uv run --with numpy --with scipy python src/SimpleReactorsBVP/analytical_reaction_diffusion_tests/generate_bvp04_reference.py

The output is a frozen numeric oracle for Rust tests; this script is never run
by the normal Rust test suite.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent / "fixtures" / "bvp04"
R = 8.314
T0 = 298.15
T_SCALE = 100.0
T_REF = 330.0
L = 0.01
RHO = 1.0
MASS_FLUX = 0.01
LAMBDA = 0.1
CP = 1000.0
DIFFUSION = 1e-4
PE_Q = MASS_FLUX * CP * L / LAMBDA
PROFILE_POINTS = 2001


@dataclass(frozen=True)
class Reaction:
    reactant: int
    stoich: tuple[float, ...]
    k_ref: float
    ea_j_mol: float


@dataclass(frozen=True)
class Scenario:
    name: str
    names: tuple[str, ...]
    molar_masses: tuple[float, ...]
    formula_atoms: tuple[dict[str, int], ...]
    inlet: tuple[float, ...]
    diffusion: tuple[float, ...]
    reactions: tuple[Reaction, ...]


SMALL_EDGES = (
    (0, 1, 0.25), (1, 2, 0.60), (2, 3, 1.00), (0, 2, 0.35), (1, 3, 0.45)
)
LARGE_EDGES = (
    (0, 1, .70), (1, 2, 1.10), (2, 3, .50), (3, 4, 1.60),
    (4, 5, .80), (5, 6, 1.30), (6, 7, .60), (7, 8, 1.90),
    (8, 9, .90), (9, 10, 1.40), (10, 11, .75), (11, 0, 1.20),
    (0, 4, .45), (1, 5, .95), (2, 6, 1.70), (3, 7, .55),
    (4, 8, 1.05), (5, 9, 1.80), (6, 10, .65), (7, 11, 1.50),
)


def chain_scenario(name: str, edges, inlet, activation_energies_kj=None) -> Scenario:
    n = len(inlet)
    reactions = []
    for j, (source, target, k_ref) in enumerate(edges):
        ea_kj = (activation_energies_kj[j] if activation_energies_kj is not None
                 else 10 + 2 * j)
        stoich = np.zeros(n)
        stoich[source] = -1.0
        stoich[target] = 1.0
        reactions.append(Reaction(source, tuple(stoich), k_ref, ea_kj * 1000.0))
    formula = {"C": 2, "H": 4, "O": 2}
    return Scenario(name, tuple(f"X{i}" for i in range(n)), (0.060052,) * n,
                    (formula,) * n, tuple(inlet), (DIFFUSION,) * n, tuple(reactions))


SCENARIOS = (
    chain_scenario("bvp04b_4_5", SMALL_EDGES, (.50, .25, .15, .10),
                   (12, 18, 24, 30, 36)),
    chain_scenario("bvp04c_12_20", LARGE_EDGES,
                   (.12, .08, .07, .09, .10, .06, .11, .05, .07, .08, .09, .08)),
    Scenario(
        "bvp04d_ozone_4_2", ("O3", "O2", "O", "N2"),
        (.047997, .031998, .015999, .028014),
        ({"O": 3}, {"O": 2}, {"O": 1}, {"N": 2}), (.35, .25, .10, .30),
        (DIFFUSION,) * 4,
        (Reaction(0, (-1, 1, 1, 0), .35, 24_000.0),
         Reaction(1, (0, -1, 2, 0), .75, 32_000.0)),
    ),
)


def temperature(z: float) -> float:
    if abs(PE_Q) < 1e-14:
        return T0
    integral = (math.exp(PE_Q * (z - 1.0)) - math.exp(-PE_Q)) / PE_Q
    q_right_internal = 1000.0 * L / T_SCALE
    theta = q_right_internal / LAMBDA * integral
    return T0 + T_SCALE * theta


def source_matrix(s: Scenario, temp: float) -> np.ndarray:
    n = len(s.names)
    w = np.zeros((n, n))
    for reaction in s.reactions:
        a = reaction.k_ref * math.exp(reaction.ea_j_mol / (R * T_REF))
        k = a * math.exp(-reaction.ea_j_mol / (R * temp))
        factor = RHO * k / s.molar_masses[reaction.reactant]
        for i, nu in enumerate(reaction.stoich):
            w[i, reaction.reactant] += s.molar_masses[i] * nu * factor
    return w


def operator(s: Scenario, z: float) -> np.ndarray:
    n = len(s.names)
    rho_d = RHO * np.asarray(s.diffusion)
    a = np.zeros((2 * n, 2 * n))
    a[:n, n:] = np.diag(1.0 / rho_d)
    a[n:, n:] = np.diag(MASS_FLUX * L / rho_d)
    a[n:, :n] = -L * L * source_matrix(s, temperature(z))
    return a


def integrate_fundamental(s: Scenario, rtol: float, max_step: float):
    n2 = 2 * len(s.names)
    eye = np.eye(n2)

    def rhs(z, flat):
        return (operator(s, z) @ flat.reshape((n2, n2))).ravel()

    return solve_ivp(
        rhs, (0.0, 1.0), eye.ravel(), method="DOP853", rtol=rtol,
        atol=rtol * 1e-3, max_step=max_step, dense_output=True,
    )


def solve_scenario(s: Scenario):
    # Establish convergence in both tolerance and maximum step before freezing.
    studies = [(1e-9, .02), (1e-11, .02), (1e-11, .01), (1e-11, .005)]
    runs = [integrate_fundamental(s, rtol, step) for rtol, step in studies]
    if not all(run.success for run in runs):
        raise RuntimeError(f"DOP853 failed for {s.name}: {[run.message for run in runs]}")
    n = len(s.names)

    def publish(run):
        n2 = 2 * n
        phi1 = run.sol(1.0).reshape((2 * n, 2 * n))
        jj = phi1[n:, n:]
        jc = phi1[n:, :n]
        initial_flux = np.linalg.solve(jj, -jc @ np.asarray(s.inlet))
        initial = np.concatenate((s.inlet, initial_flux))
        residual = np.linalg.norm(jj @ initial_flux + jc @ np.asarray(s.inlet))
        condition = float(np.linalg.cond(jj, 2))
        states = run.sol(np.linspace(0.0, 1.0, PROFILE_POINTS)).T.reshape(
            (PROFILE_POINTS, 2 * n, n2)
        )
        profiles = np.einsum("tij,j->ti", states, initial)
        return profiles, condition, residual, initial_flux

    published = [publish(run) for run in runs]
    profile_delta = float(np.max(np.abs(published[-1][0] - published[-2][0])))
    if profile_delta > 2e-8:
        raise RuntimeError(f"reference convergence failed for {s.name}: delta={profile_delta:.3e}")
    profiles, condition, boundary_residual, initial_flux = published[-1]
    n = len(s.names)
    mass_invariant = np.sum(profiles[:, n:] - MASS_FLUX * L * profiles[:, :n], axis=1)
    mass_scale = max(float(np.max(np.sum(np.abs(profiles[:, n:]), axis=1))), 1e-12)
    mass_drift = float(np.max(np.abs(mass_invariant - mass_invariant[0])) / mass_scale)
    element_drift = 0.0
    for element in {key for formula in s.formula_atoms for key in formula}:
        weights = np.asarray([
            formula.get(element, 0) / mass
            for formula, mass in zip(s.formula_atoms, s.molar_masses)
        ])
        elemental_state = profiles[:, n:] - MASS_FLUX * L * profiles[:, :n]
        elemental_invariant = elemental_state @ weights
        scale = max(float(np.max(np.sum(np.abs(elemental_state * weights), axis=1))), 1e-12)
        element_drift = max(
            element_drift,
            float(np.max(np.abs(elemental_invariant - elemental_invariant[0])) / scale),
        )
    if mass_drift > 1e-9 or element_drift > 1e-9:
        raise RuntimeError(
            f"invariant check failed for {s.name}: mass={mass_drift:.3e}, element={element_drift:.3e}"
        )
    return profiles, {
        "scenario": s.name,
        "species": list(s.names),
        "molar_masses_kg_per_mol": list(s.molar_masses),
        "formula_atom_counts": list(s.formula_atoms),
        "inlet_mass_fractions": list(s.inlet),
        "initial_internal_species_fluxes": initial_flux.tolist(),
        "condition_2_phi_jj": condition,
        "mixed_boundary_linear_residual_l2": boundary_residual,
        "mass_invariant_relative_drift": mass_drift,
        "max_element_invariant_relative_drift": element_drift,
        "convergence_max_state_delta_last_two": profile_delta,
        "study_rtol_max_step": studies,
        "scipy_version": scipy.__version__,
        "method": "DOP853",
        "profile_points": PROFILE_POINTS,
    }


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    metadata = {
        "model": "independent temperature-dependent linear reaction-diffusion BVP",
        "coordinate": "z=x/L in [0,1]",
        "temperature_formula": "T=T0 + q_phys_R*L/lambda*(exp(Pe_q*(z-1))-exp(-Pe_q))/Pe_q",
        "temperature_parameters": {"T0_K": T0, "T_scale_K": T_SCALE, "T_ref_K": T_REF,
                                   "q_phys_right_W_m2": 1000.0, "Pe_q": PE_Q},
        "scenarios": [],
    }
    for scenario in SCENARIOS:
        profiles, details = solve_scenario(scenario)
        n = len(scenario.names)
        path = ROOT / f"{scenario.name}.csv"
        header = ["z", *[f"C_{name}" for name in scenario.names],
                  *[f"J_internal_{name}" for name in scenario.names], "T_K", "q_phys_W_m2"]
        with path.open("w", newline="", encoding="ascii") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            for index, z in enumerate(np.linspace(0.0, 1.0, PROFILE_POINTS)):
                c = profiles[index, :n]
                j_internal = profiles[index, n:]
                q_phys = 1000.0 * math.exp(PE_Q * (z - 1.0))
                writer.writerow([f"{z:.17e}", *[f"{v:.17e}" for v in c],
                                 *[f"{v:.17e}" for v in j_internal],
                                 f"{temperature(z):.17e}", f"{q_phys:.17e}"])
        metadata["scenarios"].append(details)
        print(f"{scenario.name}: cond2(Phi_JJ)={details['condition_2_phi_jj']:.6e}; "
              f"boundary residual={details['mixed_boundary_linear_residual_l2']:.3e}; "
              f"convergence delta={details['convergence_max_state_delta_last_two']:.3e}")
    (ROOT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
