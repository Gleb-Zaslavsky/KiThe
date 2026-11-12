"""Independent reduced ZZ-1962-B scalar BVP reference generator.

Run from the repository root with:
    uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_zz1962b/generate_reduced_reference.py

This intentionally solves only the reduced equation, independently of KiThe.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp

ROOT = Path(__file__).resolve().parent / "fixtures"
ALPHA = 2.0
BETA = 7.0
DELTA = 10.0
ETA_MAXIMA = (10.0, 15.0, 20.0)
PROFILE_POINTS = 2001


def rhs(eta: np.ndarray, state: np.ndarray) -> np.ndarray:
    theta, slope = state
    denominator = 1.0 + DELTA * theta
    exponent = np.clip(-BETA / np.where(np.abs(denominator) < 1e-8,
                                         np.copysign(1e-8, denominator), denominator), -80.0, 80.0)
    source = (1.0 - theta) * np.exp(exponent)
    curvature = (slope - source) / ALPHA
    return np.vstack((slope, curvature))


def solve_branch(eta_max: float, seed_theta: float, tol: float = 2e-8,
                 initial_mesh_nodes: int = 101):
    eta = np.linspace(0.0, eta_max, initial_mesh_nodes)
    scale = max(0.5, eta_max / 10.0)
    shape = 1.0 - np.exp(-eta / scale)
    guess = np.vstack((seed_theta * shape, seed_theta / scale * np.exp(-eta / scale)))
    result = solve_bvp(
        rhs,
        lambda ya, yb: np.array([ya[0], yb[1]]),
        eta,
        guess,
        tol=tol,
        max_nodes=100000,
        verbose=0,
    )
    if not result.success:
        raise RuntimeError(
            f"solve_bvp failed eta_max={eta_max} seed={seed_theta}: {result.message}"
        )
    eta_out = np.linspace(0.0, eta_max, PROFILE_POINTS)
    theta, slope = result.sol(eta_out)
    if np.min(1.0 + DELTA * theta) <= 0.0:
        raise RuntimeError("reduced solution crossed the Arrhenius denominator singularity")
    collocation_residual = float(np.max(result.rms_residuals))
    bc_residual = np.array([result.y[0, 0], result.y[1, -1]])
    rate = (1.0 - theta) * np.exp(-BETA / (1.0 + DELTA * theta))
    return result, eta_out, theta, slope, rate, collocation_residual, bc_residual


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    summaries = []
    for eta_max in ETA_MAXIMA:
        branches = []
        failed_seeds = []
        for seed_theta in (0.01, 0.1, 0.3, 0.6, 0.85):
            try:
                result, eta, theta, slope, rate, ode_res, bc_res = solve_branch(
                    eta_max, seed_theta
                )
            except RuntimeError as error:
                failed_seeds.append({"seed_theta": seed_theta, "error": str(error)})
                print(f"eta_max={eta_max:g} seed={seed_theta:g} failed: {error}")
                continue
            if any(np.max(np.abs(theta - item["theta"])) < 1e-6 for item in branches):
                continue
            profile_name = f"zz1962b_eta{int(eta_max)}_branch{len(branches) + 1}.csv"
            with (ROOT / profile_name).open("w", newline="", encoding="ascii") as stream:
                writer = csv.writer(stream)
                writer.writerow(["eta", "theta", "theta_prime", "T_K", "Y_A", "Y_B", "rate_dimensionless"])
                for row in zip(eta, theta, slope):
                    eta_i, theta_i, slope_i = row
                    writer.writerow([
                        f"{eta_i:.17e}", f"{theta_i:.17e}", f"{slope_i:.17e}",
                        f"{300.0 + 3000.0 * theta_i:.17e}",
                        f"{1.0 - theta_i:.17e}", f"{theta_i:.17e}",
                        f"{(1.0 - theta_i) * np.exp(-BETA / (1.0 + DELTA * theta_i)):.17e}",
                    ])
            convergence_study = []
            for study_tol, study_nodes in ((1e-5, 51), (1e-7, 101), (2e-8, 201)):
                try:
                    study = solve_branch(eta_max, seed_theta, study_tol, study_nodes)
                    study_theta, study_slope = study[2], study[3]
                    delta = float(max(
                        np.max(np.abs(study_theta - theta)),
                        np.max(np.abs(study_slope - slope)),
                    ))
                    convergence_study.append({
                        "rtol": study_tol,
                        "initial_mesh_nodes": study_nodes,
                        "adaptive_mesh_nodes": int(study[0].x.size),
                        "max_delta_vs_frozen_profile": delta,
                        "max_collocation_rms_residual": study[5],
                    })
                except RuntimeError as error:
                    convergence_study.append({
                        "rtol": study_tol,
                        "initial_mesh_nodes": study_nodes,
                        "error": str(error),
                    })
            successful_studies = [
                item for item in convergence_study
                if "max_delta_vs_frozen_profile" in item
            ]
            if len(successful_studies) >= 2:
                fine_delta = successful_studies[-1]["max_delta_vs_frozen_profile"]
                if fine_delta > 2e-5:
                    raise RuntimeError(
                        f"reduced BVP refinement did not converge for eta_max={eta_max:g} "
                        f"seed={seed_theta}: delta={fine_delta:.3e}"
                    )
            branches.append({
                "theta": theta,
                "seed_theta": seed_theta,
                "file": profile_name,
                "mesh_nodes": int(result.x.size),
                "iterations": int(result.niter),
                "max_collocation_rms_residual": ode_res,
                "convergence_study": convergence_study,
                "bc_residuals": bc_res.tolist(),
                "theta_out": float(theta[-1]),
                "theta_min": float(np.min(theta)),
                "theta_max": float(np.max(theta)),
                "T_out_K": float(300.0 + 3000.0 * theta[-1]),
                "Y_A_out": float(1.0 - theta[-1]),
                "theta_prime_in": float(slope[0]),
                "theta_prime_out": float(slope[-1]),
                "max_rate": float(np.max(rate)),
                "eta_at_max_rate": float(eta[int(np.argmax(rate))]),
            })
        summaries.append({
            "eta_max": eta_max,
            "failed_seeds": failed_seeds,
            "branches": [{key: value for key, value in item.items() if key != "theta"}
                         for item in branches],
        })
        for branch in summaries[-1]["branches"]:
            print(
                f"eta_max={eta_max:g} branch={branch['file']} nodes={branch['mesh_nodes']} "
                f"Y_A_out={branch['Y_A_out']:.9e} T_out={branch['T_out_K']:.9f} K "
                f"eps_inlet={abs(branch['theta_prime_in']) * 2.0:.9e} "
                f"max_rate={branch['max_rate']:.9e} at eta={branch['eta_at_max_rate']:.6f} "
                f"collocation_rms={branch['max_collocation_rms_residual']:.3e}"
            )
    metadata = {
        "reference": "independent scalar reduced BVP; SciPy solve_bvp",
        "scipy_version": scipy.__version__,
        "equation": "theta_prime = alpha*theta_second + (1-theta)*exp(-beta/(1+delta*theta))",
        "first_order_rhs": "theta_second=(theta_prime-(1-theta)*exp(-beta/(1+delta*theta)))/alpha",
        "boundary_conditions": ["theta(0)=0", "theta_prime(eta_max)=0"],
        "parameters": {"alpha": ALPHA, "beta": BETA, "delta": DELTA},
        "solver_tolerance": 2e-8,
        "profile_points": PROFILE_POINTS,
        "summaries": summaries,
    }
    (ROOT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="ascii")


if __name__ == "__main__":
    main()
