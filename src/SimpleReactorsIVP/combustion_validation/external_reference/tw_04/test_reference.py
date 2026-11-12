import numpy as np
from scipy.linalg import expm

from .reference import (
    ACTIVATION_ENERGIES,
    ENTHALPIES,
    FROM,
    INITIAL_Y,
    INITIAL_QHAT,
    INITIAL_THETA,
    K_REFERENCE,
    LAMBDA,
    LENGTH,
    MASS_FLUX,
    MOLAR_MASS,
    PE_Q,
    RHO,
    REACTION_HEATS,
    SCENARIOS,
    T_SCALE,
    T_OFFSET,
    CP,
    stoichiometric_matrix,
    energy_invariant,
    reaction_heats,
    rate_constants,
    solve,
    validate_solution,
)


def test_mechanism_and_hess_law_are_closed_and_mass_balanced():
    matrix = stoichiometric_matrix()
    assert matrix.shape == (12, 20)
    assert np.allclose(np.sum(matrix, axis=0), 0.0, atol=0.0)
    assert np.allclose(ENTHALPIES @ matrix + REACTION_HEATS, 0.0, atol=0.0)
    assert ACTIVATION_ENERGIES.shape == (20,)
    assert np.isclose(np.sum(INITIAL_Y), 1.0)
    assert np.allclose(np.sum(MOLAR_MASS * matrix, axis=0), 0.0, atol=1e-15)


def test_heat_effects_scale_linearly_and_common_enthalpy_shift_cancels():
    assert np.allclose(reaction_heats(10.0), 10.0 * reaction_heats(1.0), rtol=0.0, atol=0.0)
    shifted = ENTHALPIES + 9876.5
    matrix = stoichiometric_matrix()
    assert np.allclose(-shifted @ matrix, REACTION_HEATS, rtol=0.0, atol=1e-11)


def test_all_frozen_gamma_scenarios_preserve_physical_invariants():
    for name, gamma in SCENARIOS:
        solution = solve(name, gamma, method="DOP853", rtol=1e-11, max_step=1e-3)
        diagnostics = validate_solution(solution)
        assert diagnostics["minimum_mass_fraction"] > 0.0
        assert diagnostics["mass_balance_error"] < 1e-10
        assert diagnostics["normalized_energy_invariant"] < 1e-9
        assert np.max(np.abs(energy_invariant(solution.state, gamma))) < 1e-7


def test_selected_profiles_agree_between_independent_scipy_methods():
    grid = np.linspace(0.0, 1.0, 1001)
    for name, gamma in SCENARIOS:
        explicit = solve(name, gamma, method="DOP853", rtol=1e-11, max_step=1e-3, z=grid)
        implicit = solve(name, gamma, method="Radau", rtol=1e-11, max_step=1e-3, z=grid)
        assert np.max(np.abs(explicit.state[2:] - implicit.state[2:])) < 2e-8
        assert np.max(np.abs(explicit.state[0] - implicit.state[0])) < 2e-9
        assert np.max(np.abs(explicit.state[1] - implicit.state[1])) < 2e-7


def test_zero_activation_energy_limit_matches_constant_coefficient_solution():
    gamma = 100_000.0
    zero_energies = np.zeros_like(ACTIVATION_ENERGIES)
    assert np.allclose(rate_constants(np.array([250.0, 400.0, 900.0]), zero_energies), K_REFERENCE)

    stoich = stoichiometric_matrix()
    rate_map = np.zeros((20, 12))
    for reaction, reactant in enumerate(FROM):
        rate_map[reaction, reactant] = K_REFERENCE[reaction] * RHO / MOLAR_MASS
    system = np.zeros((14, 14))
    system[0, 1] = 1.0 / LAMBDA
    system[1, 1] = PE_Q
    system[1, 2:] = -LENGTH**2 / T_SCALE * reaction_heats(gamma) @ rate_map
    system[2:, 2:] = LENGTH / MASS_FLUX * MOLAR_MASS * stoich @ rate_map
    initial = np.concatenate(([INITIAL_THETA, INITIAL_QHAT], INITIAL_Y))

    grid = np.linspace(0.0, 1.0, 101)
    actual = solve(
        "zero_activation_energy",
        gamma,
        z=grid,
        activation_energies=zero_energies,
    )
    expected = np.column_stack([expm(z * system) @ initial for z in grid])
    assert np.max(np.abs(actual.state - expected)) < 2e-9


def test_selected_gamma_values_show_measured_weak_moderate_and_strong_feedback():
    grid = np.linspace(0.0, 1.0, 1001)
    profiles = {
        name: solve(name, gamma, rtol=1e-10, z=grid, max_step=1e-3)
        for name, gamma in SCENARIOS
    }
    zero = profiles["tw_04a"].state
    weak = profiles["tw_04b"].state
    moderate = profiles["tw_04c"].state
    strong = profiles["tw_04d"].state
    d_weak = np.max(np.abs(T_OFFSET + T_SCALE * weak[0] - (T_OFFSET + T_SCALE * zero[0])))
    d_moderate = np.max(np.abs(T_OFFSET + T_SCALE * moderate[0] - (T_OFFSET + T_SCALE * zero[0])))
    d_strong = np.max(np.abs(T_OFFSET + T_SCALE * strong[0] - (T_OFFSET + T_SCALE * zero[0])))
    strong_departure = np.max(
        np.abs(
            T_SCALE * (strong[0] - zero[0])
            - (500_000.0 / 1_000.0) * T_SCALE * (weak[0] - zero[0])
        )
    )
    assert d_weak > 1e-3
    assert d_moderate > 1.0
    assert d_strong > 10.0
    assert strong_departure > 10.0
