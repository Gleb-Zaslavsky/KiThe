"""Independent numerical TW-04 model using NumPy and SciPy only."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

FLOAT = NDArray[np.float64]

R_GAS = 8.314
RHO = 1000.0
CP = 1000.0
LAMBDA = 10_000.0
LENGTH = 0.01
MASS_FLUX = 1.0
T_OFFSET = 300.0
T_SCALE = 100.0
MOLAR_MASS = 0.060052
T_REFERENCE = 400.0
INITIAL_THETA = (400.0 - T_OFFSET) / T_SCALE
INITIAL_QHAT = 5_000.0
PE_Q = MASS_FLUX * CP * LENGTH / LAMBDA
INITIAL_Y = np.array(
    [0.12, 0.08, 0.07, 0.09, 0.10, 0.06, 0.11, 0.05, 0.07, 0.08, 0.09, 0.08],
    dtype=np.float64,
)
ENTHALPIES = np.array(
    [0.0, -40.0, -25.0, -65.0, -15.0, -80.0, -35.0, -100.0, -50.0, -90.0, -60.0, -115.0],
    dtype=np.float64,
)

# from, to, k_ref [1/s], Ea [J/mol]
REACTIONS = np.array(
    [
        (0, 1, 0.70, 20_000),
        (1, 2, 1.10, 45_000),
        (2, 3, 0.50, 30_000),
        (3, 4, 1.60, 65_000),
        (4, 5, 0.80, 40_000),
        (5, 6, 1.30, 55_000),
        (6, 7, 0.60, 25_000),
        (7, 8, 1.90, 70_000),
        (8, 9, 0.90, 35_000),
        (9, 10, 1.40, 60_000),
        (10, 11, 0.75, 25_000),
        (11, 0, 1.20, 70_000),
        (0, 4, 0.45, 55_000),
        (1, 5, 0.95, 35_000),
        (2, 6, 1.70, 75_000),
        (3, 7, 0.55, 20_000),
        (4, 8, 1.05, 60_000),
        (5, 9, 1.80, 45_000),
        (6, 10, 0.65, 80_000),
        (7, 11, 1.50, 30_000),
    ],
    dtype=np.float64,
)

FROM = REACTIONS[:, 0].astype(np.int64)
TO = REACTIONS[:, 1].astype(np.int64)
K_REFERENCE = REACTIONS[:, 2]
ACTIVATION_ENERGIES = REACTIONS[:, 3]
PRE_EXPONENTIALS = K_REFERENCE * np.exp(ACTIVATION_ENERGIES / (R_GAS * T_REFERENCE))
REACTION_HEATS = ENTHALPIES[FROM] - ENTHALPIES[TO]

SCENARIOS: tuple[tuple[str, float], ...] = (
    ("tw_04a", 0.0),
    ("tw_04b", 1_000.0),
    ("tw_04c", 100_000.0),
    ("tw_04d", 500_000.0),
)


@dataclass(frozen=True)
class ReferenceSolution:
    scenario: str
    gamma: float
    z: FLOAT
    state: FLOAT
    derivative: FLOAT
    q_volume: FLOAT
    method: str
    rtol: float
    atol: float
    nfev: int


def rate_constants(
    temperature_k: FLOAT | float,
    activation_energies: FLOAT = ACTIVATION_ENERGIES,
) -> FLOAT:
    temperature = np.asarray(temperature_k, dtype=np.float64)
    pre_exponentials = K_REFERENCE * np.exp(activation_energies / (R_GAS * T_REFERENCE))
    return pre_exponentials * np.exp(-activation_energies / (R_GAS * temperature[..., np.newaxis]))


def reaction_rates(
    temperature_k: float,
    mass_fractions: FLOAT,
    activation_energies: FLOAT = ACTIVATION_ENERGIES,
) -> FLOAT:
    """Return 20 molar rates from independently computed concentrations."""
    concentrations = RHO * mass_fractions / MOLAR_MASS
    return rate_constants(temperature_k, activation_energies) * concentrations[FROM]


def stoichiometric_matrix() -> FLOAT:
    matrix = np.zeros((12, 20), dtype=np.float64)
    matrix[FROM, np.arange(20)] = -1.0
    matrix[TO, np.arange(20)] = 1.0
    return matrix


def reaction_heats(gamma: float) -> FLOAT:
    return gamma * REACTION_HEATS


def rhs(
    _z: float,
    state: FLOAT,
    gamma: float,
    activation_energies: FLOAT = ACTIVATION_ENERGIES,
) -> FLOAT:
    theta, qhat = state[:2]
    mass_fractions = state[2:]
    temperature_k = T_OFFSET + T_SCALE * theta
    rates = reaction_rates(float(temperature_k), mass_fractions, activation_energies)
    species_mass_sources = MOLAR_MASS * (stoichiometric_matrix() @ rates)
    d_y_d_z = LENGTH / MASS_FLUX * species_mass_sources
    q_volume = float(np.dot(reaction_heats(gamma), rates))
    d_theta_d_z = qhat / LAMBDA
    d_qhat_d_z = PE_Q * qhat - LENGTH**2 / T_SCALE * q_volume
    return np.concatenate(([d_theta_d_z, d_qhat_d_z], d_y_d_z))


def solve(
    scenario: str,
    gamma: float,
    *,
    method: str = "DOP853",
    rtol: float = 1e-11,
    atol: float | FLOAT = 1e-13,
    max_step: float = 1e-3,
    z: FLOAT | None = None,
    activation_energies: FLOAT = ACTIVATION_ENERGIES,
) -> ReferenceSolution:
    evaluation_grid = np.linspace(0.0, 1.0, 2001) if z is None else np.asarray(z)
    initial_state = np.concatenate(([INITIAL_THETA, INITIAL_QHAT], INITIAL_Y))
    solution = solve_ivp(
        lambda coordinate, state: rhs(coordinate, state, gamma, activation_energies),
        (0.0, 1.0),
        initial_state,
        method=method,
        t_eval=evaluation_grid,
        dense_output=True,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(f"{scenario} {method} failed: {solution.message}")
    state = np.asarray(solution.y, dtype=np.float64)
    derivative = np.column_stack(
        [
            rhs(float(coordinate), state[:, index], gamma, activation_energies)
            for index, coordinate in enumerate(solution.t)
        ]
    )
    q_volume = np.array(
        [
            np.dot(
                reaction_heats(gamma),
                reaction_rates(
                    T_OFFSET + T_SCALE * state[0, index],
                    state[2:, index],
                    activation_energies,
                ),
            )
            for index in range(state.shape[1])
        ],
        dtype=np.float64,
    )
    return ReferenceSolution(
        scenario=scenario,
        gamma=gamma,
        z=np.asarray(solution.t, dtype=np.float64),
        state=state,
        derivative=derivative,
        q_volume=q_volume,
        method=method,
        rtol=rtol,
        atol=float(np.max(np.asarray(atol))),
        nfev=solution.nfev,
    )


def energy_invariant(state: FLOAT, gamma: float) -> FLOAT:
    theta = state[0]
    qhat = state[1]
    mass_fractions = state[2:]
    g = ENTHALPIES / MOLAR_MASS
    initial_mass_fractions = INITIAL_Y.reshape(
        (INITIAL_Y.size,) + (1,) * (mass_fractions.ndim - 1)
    )
    return (
        qhat
        - INITIAL_QHAT
        - gamma
        * MASS_FLUX
        * LENGTH
        / T_SCALE
        * np.tensordot(g, mass_fractions - initial_mass_fractions, axes=(0, 0))
        - LAMBDA * PE_Q * (theta - INITIAL_THETA)
    )


def physical_temperature(state: FLOAT) -> FLOAT:
    return T_OFFSET + T_SCALE * state[0]


def physical_heat_flux(state: FLOAT) -> FLOAT:
    return state[1] * T_SCALE / LENGTH


def validate_solution(solution: ReferenceSolution) -> dict[str, float]:
    state = solution.state
    temperature = physical_temperature(state)
    mass_error = float(np.max(np.abs(np.sum(state[2:], axis=0) - 1.0)))
    minimum_fraction = float(np.min(state[2:]))
    invariant = energy_invariant(state, solution.gamma)
    initial_mass_fractions = INITIAL_Y.reshape(
        (INITIAL_Y.size,) + (1,) * (state[2:].ndim - 1)
    )
    invariant_scale = (
        np.abs(state[1])
        + abs(INITIAL_QHAT)
        + np.abs(solution.gamma * MASS_FLUX * LENGTH / T_SCALE)
        * np.sum(
            np.abs((ENTHALPIES / MOLAR_MASS).reshape((12,) + (1,) * (state[2:].ndim - 1))
                   * (state[2:] - initial_mass_fractions)),
            axis=0,
        )
        + np.abs(LAMBDA * PE_Q * (state[0] - INITIAL_THETA))
        + 1.0
    )
    normalized_energy_error = float(np.max(np.abs(invariant) / invariant_scale))
    if not np.all(np.isfinite(state)) or not np.all(np.isfinite(temperature)):
        raise AssertionError(f"{solution.scenario} contains non-finite state values")
    if float(np.min(temperature)) <= 0.0:
        raise AssertionError(f"{solution.scenario} has non-positive absolute temperature")
    if minimum_fraction < -1e-10:
        raise AssertionError(f"{solution.scenario} has negative mass fraction {minimum_fraction}")
    if mass_error > 1e-9:
        raise AssertionError(f"{solution.scenario} mass-fraction drift {mass_error}")
    if normalized_energy_error > 1e-8:
        raise AssertionError(f"{solution.scenario} energy invariant error {normalized_energy_error}")
    return {
        "temperature_min_k": float(np.min(temperature)),
        "temperature_max_k": float(np.max(temperature)),
        "minimum_mass_fraction": minimum_fraction,
        "mass_balance_error": mass_error,
        "max_energy_invariant": float(np.max(np.abs(invariant))),
        "normalized_energy_invariant": normalized_energy_error,
    }
