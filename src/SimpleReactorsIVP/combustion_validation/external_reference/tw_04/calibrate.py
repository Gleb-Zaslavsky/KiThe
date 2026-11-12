"""Explicit gamma calibration and SciPy method/refinement report."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import scipy

from .reference import SCENARIOS, physical_temperature, solve, validate_solution

OUTPUT = Path(__file__).with_name("calibration.json")
PROBES = (0.0, 1.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0, 200_000.0, 300_000.0, 500_000.0)
GRID = np.linspace(0.0, 1.0, 2001)


def profile_metrics(gamma: float, base: np.ndarray, unit: np.ndarray, current: np.ndarray) -> dict[str, float]:
    temperature = physical_temperature(current)
    base_temperature = physical_temperature(base)
    unit_temperature = physical_temperature(unit)
    return {
        "max_temperature_change_k": float(np.max(np.abs(temperature - base_temperature))),
        "max_mass_fraction_change": float(np.max(np.abs(current[2:] - base[2:]))),
        "nonlinear_temperature_departure_k": float(
            np.max(np.abs(temperature - base_temperature - gamma * (unit_temperature - base_temperature)))
        ),
    }


def main() -> None:
    reference = solve("probe_gamma_0", 0.0, z=GRID, max_step=1e-3)
    unit = solve("probe_gamma_1", 1.0, z=GRID, max_step=1e-3)
    results = []
    for gamma in PROBES:
        solution = solve(f"probe_gamma_{gamma:g}", gamma, z=GRID, max_step=1e-3)
        results.append(
            {
                "gamma": gamma,
                **profile_metrics(gamma, reference.state, unit.state, solution.state),
                **validate_solution(solution),
                "nfev": solution.nfev,
            }
        )

    selected = []
    for name, gamma in SCENARIOS:
        by_method = {}
        for method in ("DOP853", "Radau"):
            by_precision = {}
            solutions = {}
            for rtol in (1e-9, 1e-10, 1e-11):
                solution = solve(
                    name,
                    gamma,
                    method=method,
                    rtol=rtol,
                    atol=np.array([1e-11, 1e-8, *([1e-11] * 12)]),
                    max_step=np.inf,
                    z=GRID,
                )
                solutions[rtol] = solution
                by_precision[f"{rtol:.0e}"] = {
                    "nfev": solution.nfev,
                    **validate_solution(solution),
                }
            coarse = solutions[1e-10]
            fine = solutions[1e-11]
            by_method[method] = {
                "levels": by_precision,
                "rtol_1e-10_to_1e-11": {
                    "max_mass_fraction": float(np.max(np.abs(coarse.state[2:] - fine.state[2:]))),
                    "max_temperature_k": float(
                        np.max(np.abs(physical_temperature(coarse.state) - physical_temperature(fine.state)))
                    ),
                    "max_qhat": float(np.max(np.abs(coarse.state[1] - fine.state[1]))),
                },
            }
            refinement = by_method[method]["rtol_1e-10_to_1e-11"]
            if (
                refinement["max_mass_fraction"] > 3e-9
                or refinement["max_temperature_k"] > 5e-9
                or refinement["max_qhat"] > 2e-5
            ):
                raise RuntimeError(f"{name} {method}: tolerance refinement exceeds budget: {refinement}")
            selected.append((name, gamma, method, fine))

        dop = next(item[3] for item in selected if item[0] == name and item[2] == "DOP853")
        rad = next(item[3] for item in selected if item[0] == name and item[2] == "Radau")
        max_y_delta = float(np.max(np.abs(dop.state[2:] - rad.state[2:])))
        max_temperature_delta = float(
            np.max(np.abs(physical_temperature(dop.state) - physical_temperature(rad.state)))
        )
        max_qhat_delta = float(np.max(np.abs(dop.state[1] - rad.state[1])))
        if max_y_delta > 3e-9 or max_temperature_delta > 3e-9 or max_qhat_delta > 2e-5:
            raise RuntimeError(
                f"{name}: DOP853/Radau cross-check exceeds frozen-reference budget: "
                f"dY={max_y_delta:.3e}, dT={max_temperature_delta:.3e} K, "
                f"dqhat={max_qhat_delta:.3e}"
            )
        results.append(
            {
                "scenario": name,
                "gamma": gamma,
                "precision": by_method,
                "dop853_vs_radau_max_y": max_y_delta,
                "dop853_vs_radau_max_temperature_k": max_temperature_delta,
                "dop853_vs_radau_max_qhat": max_qhat_delta,
            }
        )

    payload = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "gamma_probes": results[: len(PROBES)],
        "selected_scenario_convergence": results[len(PROBES) :],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    for row in results[: len(PROBES)]:
        print(
            "gamma={gamma:g} dT={max_temperature_change_k:.6e} K "
            "dY={max_mass_fraction_change:.6e} N_T={nonlinear_temperature_departure_k:.6e} K "
            "T=[{temperature_min_k:.6f},{temperature_max_k:.6f}] K "
            "minY={minimum_mass_fraction:.6e}".format(**row)
        )
    for row in results[len(PROBES) :]:
        print(
            f"{row['scenario']} gamma={row['gamma']:g} "
            f"DOP853/Radau dT={row['dop853_vs_radau_max_temperature_k']:.3e} K "
            f"dY={row['dop853_vs_radau_max_y']:.3e}"
        )


if __name__ == "__main__":
    main()
