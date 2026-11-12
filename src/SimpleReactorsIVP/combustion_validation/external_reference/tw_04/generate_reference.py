"""Explicitly generate the checked-in TW-04 frozen profile and provenance."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
import scipy

from .reference import (
    ACTIVATION_ENERGIES,
    ENTHALPIES,
    FROM,
    K_REFERENCE,
    LENGTH,
    LAMBDA,
    MASS_FLUX,
    MOLAR_MASS,
    PRE_EXPONENTIALS,
    REACTION_HEATS,
    RHO,
    SCENARIOS,
    T_OFFSET,
    T_SCALE,
    TO,
    INITIAL_QHAT,
    INITIAL_THETA,
    INITIAL_Y,
    CP,
    solve,
    validate_solution,
)

HERE = Path(__file__).resolve().parent
FIXTURE = HERE.parent.parent / "fixtures" / "tw_04_profiles.csv"
METADATA = HERE.parent.parent / "fixtures" / "tw_04_reference_metadata.json"


def scenario_grid(gamma: float) -> tuple[np.ndarray, float]:
    base = np.linspace(0.0, 1.0, 1001)
    probe = solve("grid_probe", gamma, z=np.linspace(0.0, 1.0, 5001), max_step=2e-4)
    peak_z = float(probe.z[int(np.argmax(np.abs(probe.q_volume)))])
    local = np.array([peak_z + offset for offset in (-0.002, -0.001, 0.0, 0.001, 0.002)])
    local = local[(local >= 0.0) & (local <= 1.0)]
    return np.unique(np.concatenate((base, local))), peak_z


def main() -> None:
    rows = [
        "scenario,gamma,z,theta,qhat," + ",".join(f"Y{i}" for i in range(12))
        + "," + ",".join(f"dstate{i}" for i in range(14))
        + ",q_volume"
    ]
    metadata = {
        "specification": "TW-04",
        "reference_method": "DOP853",
        "cross_check_method": "Radau",
        "rtol": 1e-12,
        "atol": [1e-14, 1e-11, *([1e-14] * 12)],
        "max_step": 2e-4,
        "convergence_report": "../external_reference/tw_04/calibration.json",
        "convergence_rtol_levels": [1e-9, 1e-10, 1e-11],
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "physical_parameters": {
            "rho_kg_m3": RHO,
            "cp_j_kg_k": CP,
            "lambda_w_m_k": LAMBDA,
            "length_m": LENGTH,
            "mass_flux_kg_m2_s": MASS_FLUX,
            "temperature_offset_k": T_OFFSET,
            "temperature_scale_k": T_SCALE,
            "molar_mass_kg_mol": MOLAR_MASS,
            "reference_temperature_k": 400.0,
        },
        "initial_state": {
            "theta": INITIAL_THETA,
            "qhat": INITIAL_QHAT,
            "mass_fractions": INITIAL_Y.tolist(),
        },
        "species_formula": {"C": 2, "H": 4, "O": 2},
        "enthalpies_j_mol": ENTHALPIES.tolist(),
        "reaction_heats_at_gamma_1_j_mol": REACTION_HEATS.tolist(),
        "mechanism": [
            {
                "id": f"R{index + 1:02}",
                "from": int(FROM[index]),
                "to": int(TO[index]),
                "k_ref_s_inv": float(K_REFERENCE[index]),
                "activation_energy_j_mol": float(ACTIVATION_ENERGIES[index]),
                "pre_exponential_s_inv": float(PRE_EXPONENTIALS[index]),
            }
            for index in range(FROM.size)
        ],
        "scenarios": [],
    }
    for name, gamma in SCENARIOS:
        grid, peak_z = scenario_grid(gamma)
        reference = solve(
            name,
            gamma,
            method="DOP853",
            rtol=1e-12,
            atol=np.array([1e-14, 1e-11, *([1e-14] * 12)]),
            max_step=2e-4,
            z=grid,
        )
        cross = solve(
            name,
            gamma,
            method="Radau",
            rtol=1e-12,
            atol=np.array([1e-14, 1e-11, *([1e-14] * 12)]),
            max_step=2e-4,
            z=grid,
        )
        method_delta = {
            "max_mass_fraction": float(np.max(np.abs(reference.state[2:] - cross.state[2:]))),
            "max_theta": float(np.max(np.abs(reference.state[0] - cross.state[0]))),
            "max_qhat": float(np.max(np.abs(reference.state[1] - cross.state[1]))),
        }
        if max(method_delta.values()) > 2e-9:
            raise RuntimeError(f"{name}: DOP853/Radau reference mismatch: {method_delta}")
        diagnostics = validate_solution(reference)
        metadata["scenarios"].append(
            {
                "name": name,
                "gamma": gamma,
                "grid_points": int(grid.size),
                "peak_abs_q_volume_z": peak_z,
                "dop853_nfev": reference.nfev,
                "radau_nfev": cross.nfev,
                "cross_method_delta": method_delta,
                "diagnostics": diagnostics,
            }
        )
        for index, coordinate in enumerate(reference.z):
            values = [name, f"{gamma:.17e}", f"{coordinate:.17e}"]
            values.extend(f"{value:.17e}" for value in reference.state[:, index])
            values.extend(f"{value:.17e}" for value in reference.derivative[:, index])
            values.append(f"{reference.q_volume[index]:.17e}")
            rows.append(",".join(values))
        print(
            f"{name} gamma={gamma:g} points={grid.size} qvol_peak_z={peak_z:.6f} "
            f"T={diagnostics['temperature_min_k']:.6f}..{diagnostics['temperature_max_k']:.6f} K "
            f"minY={diagnostics['minimum_mass_fraction']:.6e} "
            f"energy={diagnostics['normalized_energy_invariant']:.3e}"
        )

    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text("\n".join(rows) + "\n", encoding="utf-8")
    METADATA.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {FIXTURE}")
    print(f"Wrote {METADATA}")


if __name__ == "__main__":
    main()
