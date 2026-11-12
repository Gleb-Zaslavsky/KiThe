"""Standalone SciPy solve_bvp references for the six-state MF system."""

import csv
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp


HERE = Path(__file__).resolve().parent
FIXTURES = HERE.parent / "fixtures"
BETA = 0.05
GAMMA = 0.01
LAMBDA = 0.1
DIFFUSIVITY = 1e-4


def run(
    seed_name: str,
    output_name: str,
    delta: float,
    pe: float,
    max_nodes: int,
    label: str,
):
    with (FIXTURES / seed_name).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    x = np.array([float(row["z"]) for row in rows])
    seed = np.array(
        [
            [float(row[name]) for row in rows]
            for name in ("theta", "q", "C0", "J0", "C1", "J1")
        ]
    )

    def rhs(_xi, y):
        denominator = 1.0 + BETA * y[0]
        source = delta * np.exp(y[0] / denominator) * y[2] / denominator
        return np.vstack(
            (
                y[1] / LAMBDA,
                pe * y[1] - LAMBDA * pe**2 * source,
                y[3] / DIFFUSIVITY,
                pe * y[3] + DIFFUSIVITY * GAMMA * pe**2 * source,
                y[5] / DIFFUSIVITY,
                pe * y[5] - DIFFUSIVITY * GAMMA * pe**2 * source,
            )
        )

    def boundary(left, right):
        return np.array(
            (left[0], right[1], left[2] - 1.0, right[3], left[4], right[5])
        )

    solution = solve_bvp(
        rhs, boundary, x, seed, tol=2e-7, max_nodes=max_nodes, verbose=0
    )
    if not solution.success:
        raise RuntimeError(solution.message)

    values = solution.sol(x)
    with (FIXTURES / output_name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("z", "theta", "q", "C0", "J0", "C1", "J1"))
        for column, z in enumerate(x):
            writer.writerow(
                tuple(f"{value:.17e}" for value in (z, *values[:, column]))
            )

    probe = np.linspace(x[0], x[-1], 5001)
    profile = solution.sol(probe)
    defect = solution.sol(probe, 1) - rhs(probe, profile)
    print(
        f"SciPy {scipy.__version__} six-state {label} nodes={solution.x.size} "
        f"iterations={solution.niter} max_rms={np.max(solution.rms_residuals):.9e} "
        f"max_ode_defect={np.max(np.abs(defect)):.9e} "
        f"theta_range=({np.min(profile[0]):.9e},{np.max(profile[0]):.9e}) "
        f"C0_range=({np.min(profile[2]):.9e},{np.max(profile[2]):.9e})"
    )


run(
    "mf1963_mf_high_six_state_seed.csv",
    "mf1963_mf_high_six_state_scipy_solution.csv",
    0.1,
    0.5427079731331537,
    20_000,
    "MF-HIGH",
)
run(
    "mf1963_mf_low_six_state_seed.csv",
    "mf1963_mf_low_six_state_scipy_solution.csv",
    1e-4,
    12017.324603855763,
    6_000,
    "MF-LOW clustered 2049-node seed",
)
try:
    run(
        "mf1963_mf_low_uniform_six_state_seed.csv",
        "mf1963_mf_low_uniform_six_state_scipy_solution.csv",
        1e-4,
        12017.324603855763,
        6_000,
        "MF-LOW uniform",
    )
except RuntimeError as error:
    if "maximum number of mesh nodes" not in str(error).lower():
        raise
    print(f"SciPy six-state MF-LOW uniform: SOLVE_FAIL: {error}")
else:
    raise AssertionError("uniform six-state MF-LOW control unexpectedly converged")
