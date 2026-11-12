"""Standalone SciPy counterpart for the shared two-state MF-HIGH fixture."""

import csv
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp


HERE = Path(__file__).resolve().parent
FIXTURE = HERE.parent / "fixtures" / "mf1963_mf_high_two_state_seed.csv"
SOLUTION_FIXTURE = HERE.parent / "fixtures" / "mf1963_mf_high_scipy_solution.csv"
DELTA = 0.1
BETA = 0.05
GAMMA = 0.01


def source(theta, delta):
    denominator = 1.0 + BETA * theta
    return delta * np.exp(theta / denominator) * (1.0 - GAMMA * theta) / denominator


def source_derivative(theta, delta):
    denominator = 1.0 + BETA * theta
    remaining = 1.0 - GAMMA * theta
    exponential = np.exp(theta / denominator)
    return delta * exponential * (
        remaining / denominator**3
        - GAMMA / denominator
        - BETA * remaining / denominator**2
    )


def run(
    seed_path: Path,
    output_path: Path,
    delta: float,
    label: str,
    max_nodes: int = 20_000,
) -> None:
    with seed_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    x = np.array([float(row["xi"]) for row in rows])
    seed = np.array(
        [
            [float(row["theta"]) for row in rows],
            [float(row["p"]) for row in rows],
        ]
    )
    def rhs(_x, y):
        return np.vstack((y[1], y[1] - source(y[0], delta)))

    def rhs_jac(_x, y):
        jac = np.zeros((2, 2, y.shape[1]))
        jac[0, 1] = 1.0
        jac[1, 0] = -source_derivative(y[0], delta)
        jac[1, 1] = 1.0
        return jac

    solution = solve_bvp(
        rhs,
        lambda left, right: np.array([left[0], right[1]]),
        x,
        seed,
        fun_jac=rhs_jac,
        tol=2e-7,
        max_nodes=max_nodes,
        verbose=0,
    )
    if not solution.success:
        raise RuntimeError(f"{label}: {solution.message}")

    values = solution.sol(x)
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("xi", "theta", "p"))
        for xi, theta, slope in zip(x, *values):
            writer.writerow(tuple(f"{value:.17e}" for value in (xi, theta, slope)))

    probe = np.linspace(x[0], x[-1], 5001)
    y = solution.sol(probe)
    ode_defect = solution.sol(probe, 1) - rhs(probe, y)
    seed_error = np.max(np.abs(values - seed), axis=1)
    print(
        f"SciPy two-state {label} delta={delta:g}: scipy={scipy.__version__} nodes={solution.x.size} "
        f"iterations={solution.niter} max_rms={np.max(solution.rms_residuals):.9e} "
        f"max_ode_defect={np.max(np.abs(ode_defect)):.9e} "
        f"max_seed_delta=(theta={seed_error[0]:.9e},p={seed_error[1]:.9e}) "
        f"bc=({y[0, 0]:.17e},{y[1, -1]:.17e}) "
        f"theta_range=({np.min(y[0]):.9e},{np.max(y[0]):.9e})"
    )


run(FIXTURE, SOLUTION_FIXTURE, DELTA, "MF-HIGH reference-seed")
run(
    HERE.parent / "fixtures" / "mf1963_mf_low_two_state_seed.csv",
    HERE.parent / "fixtures" / "mf1963_mf_low_scipy_solution.csv",
    1e-4,
    "MF-LOW reference-seed",
)
try:
    run(
        HERE.parent / "fixtures" / "mf1963_mf_low_uniform_two_state_seed.csv",
        HERE.parent / "fixtures" / "mf1963_mf_low_uniform_scipy_solution.csv",
        1e-4,
        "MF-LOW uniform-seed",
        max_nodes=6_000,
    )
except RuntimeError as error:
    print(f"SciPy two-state MF-LOW uniform-seed: SOLVE_FAIL: {error}")
