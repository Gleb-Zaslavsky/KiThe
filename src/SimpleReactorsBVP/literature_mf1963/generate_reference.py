"""Generate an independent SciPy BVP reference for MF-1963.

Run from the repository root:
    uv run --with numpy --with scipy --with matplotlib python src/SimpleReactorsBVP/literature_mf1963/generate_reference.py

Only NumPy and SciPy are used. This module deliberately does not import KiThe
or RustedSciThe.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_bvp, solve_ivp
from scipy.optimize import minimize_scalar


HERE = Path(__file__).resolve().parent
FIXTURE_DIR = HERE / "fixtures"
R_UNIVERSAL = 8.31446261815324
T0_K = 300.0
CP = 1000.0
RHO0 = 1.0
LAMBDA = 0.1
MASS_FLUX = 0.01
BETA = 0.05
GAMMA = 0.01
THETA_AD = 1.0 / GAMMA
T_SCALE = R_UNIVERSAL * T0_K**2 / (R_UNIVERSAL * T0_K / BETA)
ENERGY_ACTIVATION = R_UNIVERSAL * T0_K / BETA
HEAT_RELEASE = CP * T0_K * BETA / GAMMA
CASES = (
    ("MF-LOW", 1e-4, "long induction; weak inlet heat flux"),
    ("MF-HIGH", 0.1, "compact intense front near the inlet"),
)
TAIL_STUDIES = (1e-5, 1e-8, 1e-10)
MESH_TOLERANCES = (2e-5, 2e-7)
REFERENCE_TOLERANCE = MESH_TOLERANCES[-1]
PROFILE_BASE_POINTS = 801


@dataclass
class SolveResult:
    name: str
    delta: float
    tail_epsilon: float
    tolerance: float
    solution: object
    domain_length: float
    seed_node_count: int


def source(theta: np.ndarray | float, delta: float) -> np.ndarray | float:
    denominator = 1.0 + BETA * theta
    return (
        delta
        * np.exp(theta / denominator)
        * (1.0 - GAMMA * theta)
        / denominator
    )


def source_derivative(theta: np.ndarray | float, delta: float) -> np.ndarray | float:
    denominator = 1.0 + BETA * theta
    remaining = 1.0 - GAMMA * theta
    exponential = np.exp(theta / denominator)
    return delta * exponential * (
        remaining / denominator**3
        - GAMMA / denominator
        - BETA * remaining / denominator**2
    )


def stable_eigenvalue(delta: float) -> float:
    derivative_at_equilibrium = float(source_derivative(THETA_AD, delta))
    discriminant = 1.0 - 4.0 * derivative_at_equilibrium
    if derivative_at_equilibrium >= 0.0 or discriminant <= 0.0:
        raise ValueError("the adiabatic equilibrium has no real stable eigenvalue")
    return 0.5 * (1.0 - np.sqrt(discriminant))


def stable_manifold_seed(delta: float, tail_epsilon: float):
    eigenvalue = stable_eigenvalue(delta)
    initial = np.array(
        [THETA_AD - tail_epsilon, -eigenvalue * tail_epsilon], dtype=float
    )

    def rhs(_xi: float, state: np.ndarray) -> np.ndarray:
        theta, slope = state
        return np.array([slope, slope - source(theta, delta)])

    def inlet_event(_xi: float, state: np.ndarray) -> float:
        return float(state[0])

    inlet_event.terminal = True
    inlet_event.direction = -1.0
    trajectory = solve_ivp(
        rhs,
        (0.0, -1e9),
        initial,
        method="DOP853",
        events=inlet_event,
        rtol=2e-11,
        atol=2e-13,
        max_step=20.0,
    )
    if not trajectory.success or not trajectory.t_events[0].size:
        raise RuntimeError(
            f"stable-manifold seed did not reach theta=0 for delta={delta:g}: "
            f"{trajectory.message}"
        )

    domain_length = -float(trajectory.t[-1])
    mesh = trajectory.t[::-1] + domain_length
    guess = trajectory.y[:, ::-1].copy()
    if np.any(np.diff(mesh) <= 0.0):
        raise RuntimeError("stable-manifold mesh is not strictly increasing")
    guess[0, 0] = 0.0
    guess[0, -1] = THETA_AD
    guess[1, -1] = 0.0
    return mesh, guess, domain_length


def solve_case(name: str, delta: float, tail_epsilon: float, tolerance: float) -> SolveResult:
    mesh, guess, domain_length = stable_manifold_seed(delta, tail_epsilon)

    def rhs(_xi: np.ndarray, state: np.ndarray) -> np.ndarray:
        theta, slope = state
        return np.vstack((slope, slope - source(theta, delta)))

    def rhs_jac(_xi: np.ndarray, state: np.ndarray) -> np.ndarray:
        theta = state[0]
        result = np.zeros((2, 2, theta.size), dtype=float)
        result[0, 1, :] = 1.0
        result[1, 0, :] = -source_derivative(theta, delta)
        result[1, 1, :] = 1.0
        return result

    result = solve_bvp(
        rhs,
        lambda left, right: np.array([left[0], right[1]]),
        mesh,
        guess,
        fun_jac=rhs_jac,
        tol=tolerance,
        max_nodes=100_000,
        verbose=0,
    )
    if not result.success:
        raise RuntimeError(
            f"SciPy solve_bvp failed for {name}, tail={tail_epsilon:g}, "
            f"tol={tolerance:g}: {result.message}"
        )
    return SolveResult(
        name=name,
        delta=delta,
        tail_epsilon=tail_epsilon,
        tolerance=tolerance,
        solution=result,
        domain_length=domain_length,
        seed_node_count=int(mesh.size),
    )


def solution_diagnostics(case: SolveResult) -> dict[str, float | int]:
    sol = case.solution
    length = case.domain_length
    left_theta, left_slope = sol.sol(0.0)
    right_theta, right_slope = sol.sol(length)
    grid = np.linspace(0.0, length, 5001)
    theta_grid, slope_grid = sol.sol(grid)
    source_grid = source(theta_grid, case.delta)
    peak_index = int(np.argmax(source_grid))
    lo = grid[max(0, peak_index - 1)]
    hi = grid[min(grid.size - 1, peak_index + 1)]
    peak = minimize_scalar(
        lambda xi: -float(source(sol.sol(xi)[0], case.delta)),
        bounds=(float(lo), float(hi)),
        method="bounded",
        options={"xatol": 1e-13},
    )
    integral, integral_error = quad(
        lambda xi: float(source(sol.sol(xi)[0], case.delta)),
        0.0,
        length,
        epsabs=2e-9,
        epsrel=2e-10,
        limit=1000,
    )
    identity = float(right_theta - left_theta - right_slope + left_slope)
    xi_theta_90 = float(
        np.interp(90.0, theta_grid, grid)
    ) if theta_grid[-1] >= 90.0 else float("nan")
    convective_xi_theta_90 = float(
        quad(
            lambda theta: 1.0 / float(source(theta, case.delta)),
            0.0,
            90.0,
            epsabs=2e-8,
            epsrel=2e-10,
            limit=500,
        )[0]
    )
    dimensional = dimensional_parameters(case.delta)
    return {
        "domain_length_xi": length,
        "domain_length_m": length * float(dimensional["length_scale_m"]),
        "seed_mesh_nodes": case.seed_node_count,
        "adaptive_mesh_nodes": int(sol.x.size),
        "collocation_iterations": int(sol.niter),
        "max_rms_collocation_residual": float(np.max(sol.rms_residuals)),
        "bc_theta_inlet": float(left_theta),
        "bc_slope_outlet": float(right_slope),
        "theta_outlet": float(right_theta),
        "tail_conversion_error": float(abs(1.0 - GAMMA * right_theta)),
        "inlet_slope": float(left_slope),
        "inlet_heat_flux_W_per_m2": float(T_SCALE * MASS_FLUX * CP * left_slope),
        "xi_theta_90": xi_theta_90,
        "x_theta_90_m": xi_theta_90 * float(dimensional["length_scale_m"]),
        "convective_approx_xi_theta_90": convective_xi_theta_90,
        "convective_approx_relative_delta_at_theta_90": float(
            (xi_theta_90 - convective_xi_theta_90) / convective_xi_theta_90
        ),
        "xi_source_max": float(peak.x),
        "x_source_max_m": float(peak.x) * float(dimensional["length_scale_m"]),
        "source_max": float(-peak.fun),
        "source_integral_quad": float(integral),
        "source_integral_quad_error_estimate": float(integral_error),
        "source_integral_ode_identity": identity,
        "source_integral_identity_abs_error": float(abs(integral - identity)),
        "temperature_inlet_K": T0_K,
        "temperature_outlet_K": float(T0_K * (1.0 + BETA * right_theta)),
        "reactant_fraction_outlet": float(1.0 - GAMMA * right_theta),
        "arrhenius_k0_per_s": dimensional["k0_per_s"],
        "mass_heat_release_J_per_kg": HEAT_RELEASE,
        "thermal_diffusivity_density_product_kg_per_m_s": LAMBDA / CP,
    }


def dimensional_parameters(delta: float) -> dict[str, float]:
    k_at_reference = delta * CP * GAMMA * MASS_FLUX**2 / (LAMBDA * RHO0)
    k0 = k_at_reference * np.exp(ENERGY_ACTIVATION / (R_UNIVERSAL * T0_K))
    return {
        "T0_K": T0_K,
        "cp_J_per_kg_K": CP,
        "rho0_kg_per_m3": RHO0,
        "lambda_W_per_m_K": LAMBDA,
        "mass_flux_kg_per_m2_s": MASS_FLUX,
        "E_J_per_mol": ENERGY_ACTIVATION,
        "Q_J_per_kg": HEAT_RELEASE,
        "rhoD_kg_per_m_s": LAMBDA / CP,
        "k_at_T0_per_s": float(k_at_reference),
        "k0_per_s": float(k0),
        "temperature_scale_K": float(T_SCALE),
        "length_scale_m": float(LAMBDA / (MASS_FLUX * CP)),
    }


def profile_coordinates(case: SolveResult, diagnostics: dict[str, float | int]):
    sol = case.solution
    length = case.domain_length
    width = max(length * 1e-7, min(length * 0.02, 1.0 / max(diagnostics["source_max"], 1e-30)))
    center = float(diagnostics["xi_source_max"])
    local = np.linspace(max(0.0, center - 12.0 * width), min(length, center + 12.0 * width), 2001)
    coordinates = np.unique(np.concatenate((
        np.linspace(0.0, length, PROFILE_BASE_POINTS),
        case.solution.x,
        local,
        np.array([0.0, length]),
    )))
    return coordinates


def write_case_profile(case: SolveResult, diagnostics: dict[str, float | int]) -> str:
    coordinates = profile_coordinates(case, diagnostics)
    theta, slope = case.solution.sol(coordinates)
    rate = source(theta, case.delta)
    temperature = T0_K * (1.0 + BETA * theta)
    reactant = 1.0 - GAMMA * theta
    safe_name = case.name.lower().replace("-", "_")
    filename = f"mf1963_{safe_name}.csv"
    with (FIXTURE_DIR / filename).open("w", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream)
        writer.writerow(["xi", "theta", "dtheta_dxi", "reaction_source", "temperature_K", "Y_A"])
        for values in zip(coordinates, theta, slope, rate, temperature, reactant):
            writer.writerow([f"{float(value):.17e}" for value in values])
    return filename


def write_case_plot(case: SolveResult, diagnostics: dict[str, float | int]) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    matplotlib.rcParams["svg.hashsalt"] = "kith-mf1963-reference"
    length = case.domain_length
    xi = np.linspace(0.0, length, 4001)
    theta, _slope = case.solution.sol(xi)
    temperature = T0_K * (1.0 + BETA * theta)
    source_values = source(theta, case.delta)
    above = np.flatnonzero(source_values >= 0.1 * diagnostics["source_max"])
    if above.size:
        x_lo, x_hi = xi[above[0]], xi[above[-1]]
        margin = max(0.15 * (x_hi - x_lo), length * 1e-8)
        rate_x = np.linspace(max(0.0, x_lo - margin), min(length, x_hi + margin), 2500)
    else:
        rate_x = xi
    rate_theta = case.solution.sol(rate_x)[0]
    rate_values = source(rate_theta, case.delta)

    fig, axes = plt.subplots(2, 1, figsize=(8.2, 6.0), constrained_layout=True)
    axes[0].plot(xi, temperature, color="#b43c36", linewidth=1.5)
    axes[0].set_ylabel("Temperature T [K]")
    axes[0].set_xlabel("Dimensionless axial coordinate xi")
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(rate_x, rate_values, color="#236a8a", linewidth=1.5)
    axes[1].set_ylabel("Dimensionless source S(theta)")
    axes[1].set_xlabel("Dimensionless axial coordinate xi (front window)")
    axes[1].grid(True, alpha=0.25)
    fig.suptitle(f"{case.name}: delta={case.delta:g}, beta={BETA:g}, gamma={GAMMA:g}")
    filename = f"{case.name.lower().replace('-', '_')}_profiles.svg"
    fig.savefig(FIXTURE_DIR / filename, format="svg", metadata={"Date": None})
    plt.close(fig)
    return filename


def solve_studies(name: str, delta: float) -> tuple[SolveResult, dict, list[dict]]:
    reference = solve_case(name, delta, TAIL_STUDIES[-1], REFERENCE_TOLERANCE)
    reference_info = solution_diagnostics(reference)
    studies: list[dict] = []

    for epsilon in TAIL_STUDIES[:-1]:
        domain_case = solve_case(name, delta, epsilon, REFERENCE_TOLERANCE)
        common_end = min(reference.domain_length, domain_case.domain_length)
        common_x = np.linspace(0.0, common_end, 3001)
        ref_profile = reference.solution.sol(common_x)
        studied_profile = domain_case.solution.sol(common_x)
        max_theta_delta = float(np.max(np.abs(ref_profile[0] - studied_profile[0])))
        max_slope_delta = float(np.max(np.abs(ref_profile[1] - studied_profile[1])))
        studies.append({
            "kind": "finite_domain_tail",
            "tail_epsilon": epsilon,
            "domain_length_xi": domain_case.domain_length,
            "adaptive_mesh_nodes": int(domain_case.solution.x.size),
            "common_domain_max_theta_delta_vs_reference": max_theta_delta,
            "common_domain_max_slope_delta_vs_reference": max_slope_delta,
        })

    coarse = solve_case(name, delta, TAIL_STUDIES[-1], MESH_TOLERANCES[0])
    common_x = np.linspace(0.0, reference.domain_length, 3001)
    coarse_profile = coarse.solution.sol(common_x)
    fine_profile = reference.solution.sol(common_x)
    studies.append({
        "kind": "collocation_tolerance",
        "tolerance": MESH_TOLERANCES[0],
        "adaptive_mesh_nodes": int(coarse.solution.x.size),
        "reference_tolerance": REFERENCE_TOLERANCE,
        "reference_mesh_nodes": int(reference.solution.x.size),
        "max_theta_delta_vs_reference": float(np.max(np.abs(coarse_profile[0] - fine_profile[0]))),
        "max_slope_delta_vs_reference": float(np.max(np.abs(coarse_profile[1] - fine_profile[1]))),
    })
    return reference, reference_info, studies


def assert_numerical_contract(case: SolveResult, info: dict) -> None:
    sol = case.solution
    sample_x = np.linspace(0.0, case.domain_length, 4001)
    theta, slope = sol.sol(sample_x)
    denominator = 1.0 + BETA * theta
    reactant = 1.0 - GAMMA * theta
    assert np.isfinite(theta).all() and np.isfinite(slope).all()
    assert np.isfinite(source(theta, case.delta)).all()
    assert np.all(denominator > 0.0)
    assert np.min(reactant) >= -2e-7 and np.max(reactant) <= 1.0 + 2e-7
    assert abs(info["bc_theta_inlet"]) <= 2e-7
    assert abs(info["bc_slope_outlet"]) <= 2e-7
    assert info["max_rms_collocation_residual"] <= max(5.0 * case.tolerance, 2e-6)
    assert info["source_integral_identity_abs_error"] <= max(
        2e-6, 50.0 * info["source_integral_quad_error_estimate"]
    )
    assert info["source_max"] > 0.0


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    case_records = []
    for name, delta, regime in CASES:
        reference, info, studies = solve_studies(name, delta)
        assert_numerical_contract(reference, info)
        profile_name = write_case_profile(reference, info)
        plot_name = write_case_plot(reference, info)
        case_records.append({
            "name": name,
            "delta": delta,
            "regime": regime,
            "scipy_status": int(reference.solution.status),
            "scipy_message": reference.solution.message,
            "solver_tolerance": reference.tolerance,
            "tail_epsilon": reference.tail_epsilon,
            "profile_csv": profile_name,
            "plot_svg": plot_name,
            "profile_rows": int(len(profile_coordinates(reference, info))),
            "dimensionless_parameters": {
                "beta": BETA,
                "gamma": GAMMA,
                "reaction_order": 1,
                "theta_adiabatic": THETA_AD,
            },
            "dimensional_parameters": dimensional_parameters(delta),
            "diagnostics": info,
            "convergence_studies": studies,
        })
        print(
            f"{name} delta={delta:g} xi_max={info['domain_length_xi']:.6g} "
            f"nodes={info['adaptive_mesh_nodes']} theta_inlet_slope={info['inlet_slope']:.9g} "
            f"xi_source_max={info['xi_source_max']:.6g} "
            f"source_max={info['source_max']:.6g} "
            f"integral_error={info['source_integral_identity_abs_error']:.3e} "
            f"collocation_rms={info['max_rms_collocation_residual']:.3e}"
        )

    metadata = {
        "reference": "Independent Merzhanov-Filonenko first-order flow BVP; SciPy solve_bvp",
        "scipy_version": scipy.__version__,
        "numpy_version": np.__version__,
        "universal_gas_constant_J_per_mol_K": R_UNIVERSAL,
        "equations": [
            "theta_prime = p",
            "p_prime = p - source(theta)",
            "source = delta*exp(theta/(1+beta*theta))*(1-gamma*theta)/(1+beta*theta)",
        ],
        "boundary_conditions": ["theta(0)=0", "p(xi_max)=0"],
        "solver": {
            "method": "scipy.integrate.solve_bvp",
            "jacobian": "analytic derivative of the dimensionless RHS",
            "reference_tolerance": REFERENCE_TOLERANCE,
            "mesh_convergence_tolerances": MESH_TOLERANCES,
            "tail_epsilon_studies": TAIL_STUDIES,
        },
        "fixed_parameters": {
            "T0_K": T0_K,
            "cp_J_per_kg_K": CP,
            "rho0_kg_per_m3": RHO0,
            "lambda_W_per_m_K": LAMBDA,
            "mass_flux_kg_per_m2_s": MASS_FLUX,
            "beta": BETA,
            "gamma": GAMMA,
            "adiabatic_temperature_K": T0_K * (1.0 + BETA * THETA_AD),
            "rhoD_kg_per_m_s": LAMBDA / CP,
        },
        "cases": case_records,
    }
    (FIXTURE_DIR / "metadata.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="ascii"
    )


if __name__ == "__main__":
    main()
