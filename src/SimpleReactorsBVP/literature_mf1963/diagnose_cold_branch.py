"""Independent MF-1963 cold-branch continuation; imports neither KiThe nor RST.

Run from the repository root:
    uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/diagnose_cold_branch.py

The continuous solve uses SciPy collocation. The discrete solve is an explicit
scalar trapezoid system on the unchanged MF-HIGH interval, with an optional
pseudo-arclength corrector. Results are written next to this script.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp
from scipy.sparse import bmat, csc_matrix, diags, lil_matrix
from scipy.sparse.linalg import svds, spsolve


HERE = Path(__file__).resolve().parent
XI_MAX = 0.5427079731331537
BETA = 0.05
GAMMA = 0.01
DELTA_START = 0.1
STATE_SCALE = np.array([0.05, 0.05])


def source(theta: np.ndarray | float, delta: float) -> np.ndarray | float:
    den = 1.0 + BETA * theta
    with np.errstate(over="ignore", invalid="ignore"):
        return delta * np.exp(theta / den) * (1.0 - GAMMA * theta) / den


def source_dtheta(theta: np.ndarray | float, delta: float) -> np.ndarray | float:
    den = 1.0 + BETA * theta
    rem = 1.0 - GAMMA * theta
    return delta * np.exp(theta / den) * (
        rem / den**3 - GAMMA / den - BETA * rem / den**2
    )


def source_ddelta(theta: np.ndarray | float) -> np.ndarray | float:
    den = 1.0 + BETA * theta
    return np.exp(theta / den) * (1.0 - GAMMA * theta) / den


def cold_guess(delta: float, mesh: np.ndarray) -> np.ndarray:
    # Exact solution of the linearization S(theta)=delta at theta=0.
    p = delta * (1.0 - np.exp(mesh - XI_MAX))
    theta = np.zeros_like(mesh)
    theta[1:] = np.cumsum(0.5 * np.diff(mesh) * (p[:-1] + p[1:]))
    return np.vstack((theta, p))


def hot_guess(delta: float, mesh: np.ndarray) -> np.ndarray:
    reference = np.genfromtxt(HERE / "fixtures" / "mf1963_mf_high.csv", delimiter=",", names=True)
    theta = np.interp(mesh, reference["xi"], reference["theta"])
    p = np.interp(mesh, reference["xi"], reference["dtheta_dxi"])
    if delta != 0.1:
        raise ValueError("the frozen hot-profile seed is defined at delta=0.1")
    return np.vstack((theta, p))


def collocation_continue() -> list[dict[str, float | int | str]]:
    x = np.linspace(0.0, XI_MAX, 121)
    delta = DELTA_START
    y = cold_guess(delta, x)
    rows: list[dict[str, float | int | str]] = []
    targets = [DELTA_START, 0.2, 0.5, 1.0, 2.0, 3.0, 3.8, 4.0, 4.1,
               4.12, 4.14, 4.142]
    targets.extend(round(4.142 + i * 0.0001, 7) for i in range(1, 11))
    targets.extend([4.15, 4.2, 4.3, 4.5, 5.0, 7.5, 10.0])

    def solve_at(value: float, mesh: np.ndarray, guess: np.ndarray):
        def rhs(_x: np.ndarray, state: np.ndarray) -> np.ndarray:
            theta, p = state
            return np.vstack((p, p - source(theta, value)))

        def rhs_jac(_x: np.ndarray, state: np.ndarray) -> np.ndarray:
            result = np.zeros((2, 2, state.shape[1]))
            result[0, 1] = 1.0
            result[1, 0] = -source_dtheta(state[0], value)
            result[1, 1] = 1.0
            return result

        return solve_bvp(
            rhs,
            lambda left, right: np.array([left[0], right[1]]),
            mesh,
            guess,
            fun_jac=rhs_jac,
            bc_jac=lambda *_: (
                np.array([[1.0, 0.0], [0.0, 0.0]]),
                np.array([[0.0, 0.0], [0.0, 1.0]]),
            ),
            tol=1e-8,
            max_nodes=100_000,
        )

    previous = None
    for target in targets:
        if previous is not None:
            x = previous.x
            y = previous.y
        result = solve_at(target, x, y)
        if not result.success:
            rows.append({"delta": target, "status": "failed", "message": result.message})
            break
        previous = result
        sample_x = np.linspace(0.0, XI_MAX, 2001)
        sample = result.sol(sample_x)
        row = {
            "delta": float(target),
            "status": "ok",
            "nodes": int(result.x.size),
            "iterations": int(result.niter),
            "theta_max": float(np.max(sample[0])),
            "p_inlet": float(result.sol(0.0)[1]),
            "temperature_max_K": float(300.0 + 15.0 * np.max(sample[0])),
            "bc_inf": float(max(abs(result.sol(0.0)[0]), abs(result.sol(XI_MAX)[1]))),
            "max_rms_collocation_residual": float(np.max(result.rms_residuals)),
        }
        rows.append(row)
        print("COLLOCATION " + " ".join(f"{k}={v}" for k, v in row.items()))
    return rows


def trapezoid_residual(u: np.ndarray, delta: float, intervals: int) -> np.ndarray:
    y = u.reshape(intervals + 1, 2)
    theta, p = y[:, 0], y[:, 1]
    h = XI_MAX / intervals
    f = p - source(theta, delta)
    residual = np.empty_like(u)
    r = residual.reshape(intervals + 1, 2)
    r[0, 0] = theta[0]
    r[0, 1] = theta[1] - theta[0] - 0.5 * h * (p[0] + p[1])
    for i in range(1, intervals):
        r[i, 0] = p[i] - p[i - 1] - 0.5 * h * (f[i - 1] + f[i])
        r[i, 1] = theta[i + 1] - theta[i] - 0.5 * h * (p[i] + p[i + 1])
    r[intervals, 0] = p[intervals] - p[intervals - 1] - 0.5 * h * (f[intervals - 1] + f[intervals])
    r[intervals, 1] = p[intervals]
    return residual


def trapezoid_parameter_derivative(u: np.ndarray, intervals: int) -> np.ndarray:
    """Exact dF/ddelta in the same interleaved row order as trapezoid_residual."""
    theta = u.reshape(intervals + 1, 2)[:, 0]
    h = XI_MAX / intervals
    g = source_ddelta(theta)
    derivative = np.zeros((intervals + 1, 2), dtype=float)
    # Row 0 is the left theta boundary condition. Each p balance is stored
    # in column 0 of rows 1..N, including the final interval at row N.
    derivative[1:, 0] = 0.5 * h * (g[:-1] + g[1:])
    return derivative.reshape(-1)


def trapezoid_jacobian(u: np.ndarray, delta: float, intervals: int) -> csc_matrix:
    y = u.reshape(intervals + 1, 2)
    theta, p = y[:, 0], y[:, 1]
    h = XI_MAX / intervals
    jac = lil_matrix((u.size, u.size))
    jac[0, 0] = 1.0
    # First interval's theta equation is at row 1, with its left boundary row at 0.
    jac[1, 0] = -1.0
    jac[1, 2] = 1.0
    jac[1, 1] = -0.5 * h
    jac[1, 3] = -0.5 * h
    for i in range(1, intervals):
        row_p = 2 * i
        row_theta = row_p + 1
        left, right = 2 * (i - 1), 2 * i
        jac[row_p, left] = 0.5 * h * source_dtheta(theta[i - 1], delta)
        jac[row_p, left + 1] = -1.0 - 0.5 * h
        jac[row_p, right] = 0.5 * h * source_dtheta(theta[i], delta)
        jac[row_p, right + 1] = 1.0 - 0.5 * h
        jac[row_theta, right] = -1.0
        jac[row_theta, right + 1] = -0.5 * h
        jac[row_theta, right + 2] = 1.0
        jac[row_theta, right + 3] = -0.5 * h
    row_p = 2 * intervals
    left, right = 2 * (intervals - 1), 2 * intervals
    jac[row_p, left] = 0.5 * h * source_dtheta(theta[-2], delta)
    jac[row_p, left + 1] = -1.0 - 0.5 * h
    jac[row_p, right] = 0.5 * h * source_dtheta(theta[-1], delta)
    jac[row_p, right + 1] = 1.0 - 0.5 * h
    jac[2 * intervals + 1, right + 1] = 1.0
    return jac.tocsc()


def newton_correct(u: np.ndarray, delta: float, intervals: int, tol: float = 1e-11):
    for iteration in range(35):
        residual = trapezoid_residual(u, delta, intervals)
        norm = np.linalg.norm(residual, ord=np.inf)
        if norm < tol:
            return u, True, iteration, norm
        step = spsolve(trapezoid_jacobian(u, delta, intervals), -residual)
        alpha = 1.0
        base = np.linalg.norm(residual)
        accepted = False
        while alpha >= 2.0**-20:
            candidate = u + alpha * step
            trial = trapezoid_residual(candidate, delta, intervals)
            if np.isfinite(trial).all() and np.linalg.norm(trial) < base:
                u = candidate
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return u, False, iteration + 1, norm
    return u, False, 35, float(np.linalg.norm(trapezoid_residual(u, delta, intervals), ord=np.inf))


def branch_row(u: np.ndarray, delta: float, intervals: int, status: str, iterations: int, spectral: bool = False) -> dict[str, float | int | str]:
    residual = trapezoid_residual(u, delta, intervals)
    y = u.reshape(intervals + 1, 2)
    jac = trapezoid_jacobian(u, delta, intervals)
    row = {
        "delta": float(delta),
        "status": status,
        "intervals": intervals,
        "iterations": iterations,
        "theta_max": float(np.max(y[:, 0])),
        "temperature_max_K": float(300.0 + 15.0 * np.max(y[:, 0])),
        "p_inlet": float(y[0, 1]),
        "residual_inf": float(np.linalg.norm(residual, ord=np.inf)),
        "residual_l2": float(np.linalg.norm(residual)),
    }
    if spectral:
        try:
            sigma_min = float(svds(jac, k=1, which="SM", return_singular_vectors=False, tol=1e-7, maxiter=20_000)[0])
            sigma_max = float(svds(jac, k=1, which="LM", return_singular_vectors=False, tol=1e-7, maxiter=20_000)[0])
        except Exception:
            sigma_min = float("nan")
            sigma_max = float("nan")
        row.update({"sigma_min_raw_J": sigma_min, "sigma_max_raw_J": sigma_max,
                    "cond_raw_J_est": sigma_max / sigma_min if sigma_min > 0.0 else float("inf")})
    return row


def fixed_scaled_spectrum(u: np.ndarray, delta: float, intervals: int) -> dict[str, float]:
    """Spectrum with fixed state scales and mesh-consistent residual row scales."""
    jac = trapezoid_jacobian(u, delta, intervals)
    h = XI_MAX / intervals
    state_scales = np.tile(STATE_SCALE, intervals + 1)
    row_scales = np.full(2 * (intervals + 1), 1.0 / h)
    row_scales[0] = row_scales[-1] = 1.0
    scaled = diags(row_scales) @ jac @ diags(state_scales)

    def extremes(matrix: csc_matrix) -> tuple[float, float]:
        small = float(svds(matrix, k=1, which="SM", return_singular_vectors=False,
                           tol=1e-7, maxiter=30_000)[0])
        large = float(svds(matrix, k=1, which="LM", return_singular_vectors=False,
                           tol=1e-7, maxiter=30_000)[0])
        return small, large

    raw_min, raw_max = extremes(jac)
    scaled_min, scaled_max = extremes(scaled.tocsc())
    return {
        "delta": float(delta),
        "sigma_min_raw": raw_min,
        "sigma_max_raw": raw_max,
        "condition_raw": raw_max / raw_min,
        "sigma_min_fixed_scaled": scaled_min,
        "sigma_max_fixed_scaled": scaled_max,
        "condition_fixed_scaled": scaled_max / scaled_min,
    }


def discrete_continue(intervals: int, deltas: list[float]) -> list[dict[str, float | int | str]]:
    mesh = np.linspace(0.0, XI_MAX, intervals + 1)
    u = cold_guess(DELTA_START, mesh).T.reshape(-1)
    rows = []
    for delta in deltas:
        if delta == DELTA_START:
            seed = u
        else:
            seed = u.copy()
        root, ok, iterations, norm = newton_correct(seed, delta, intervals)
        if not ok:
            rows.append(branch_row(root, delta, intervals, "failed", iterations) | {"failure_norm_inf": norm})
            print("TRAPEZOID " + " ".join(f"{k}={v}" for k, v in rows[-1].items()))
            break
        u = root
        spectral = intervals == 301 and delta in {0.1, 4.0, 4.12, 4.142, 4.143, 4.15}
        rows.append(branch_row(root, delta, intervals, "ok", iterations, spectral=spectral))
        print("TRAPEZOID " + " ".join(f"{k}={v}" for k, v in rows[-1].items()))
    return rows


def discrete_hot_root(intervals: int, deltas: list[float]) -> list[dict[str, float | int | str]]:
    mesh = np.linspace(0.0, XI_MAX, intervals + 1)
    u = hot_guess(0.1, mesh).T.reshape(-1)
    rows = []
    for delta in deltas:
        root, ok, iterations, norm = newton_correct(u.copy(), delta, intervals)
        if not ok:
            rows.append(branch_row(root, delta, intervals, "failed", iterations) | {"failure_norm_inf": norm})
            break
        u = root
        rows.append(branch_row(root, delta, intervals, "ok", iterations))
    return rows


def pseudo_arclength(intervals: int, points: int = 1200, ds: float = 0.03):
    mesh = np.linspace(0.0, XI_MAX, intervals + 1)
    delta = 4.0
    u = cold_guess(delta, mesh).T.reshape(-1)
    u, ok, _, _ = newton_correct(u, delta, intervals)
    if not ok:
        raise RuntimeError("delta=4 cold point did not solve on the trapezoid grid")
    scales = np.tile(STATE_SCALE, intervals + 1)
    weights = np.full(u.size, 1.0 / (intervals + 1))
    rows = []
    fdelta = trapezoid_parameter_derivative(u, intervals)
    jac = trapezoid_jacobian(u, delta, intervals)
    tangent_q_delta = spsolve(jac @ csc_matrix(np.diag(scales)), -fdelta)
    tangent = np.append(tangent_q_delta, 1.0)
    tangent /= np.sqrt(np.dot(weights, tangent[:-1] ** 2) + tangent[-1] ** 2)

    for _ in range(points):
        step_ds = ds
        q = u / scales
        predictor = np.append(q, delta) + step_ds * tangent
        qk = predictor[:-1].copy()
        dk = float(predictor[-1])
        success = False
        corrector_residual_inf = float("inf")
        for _iteration in range(18):
            uk = qk * scales
            f = trapezoid_residual(uk, dk, intervals)
            fd = trapezoid_parameter_derivative(uk, intervals)
            difference = np.append(qk, dk) - predictor
            arc = float(np.dot(weights * tangent[:-1], difference[:-1]) + tangent[-1] * difference[-1])
            augmented = np.append(f, arc)
            corrector_residual_inf = float(np.linalg.norm(augmented, ord=np.inf))
            if corrector_residual_inf < 2e-10:
                success = True
                break
            j_scaled = trapezoid_jacobian(uk, dk, intervals) @ csc_matrix(np.diag(scales))
            border = bmat([[j_scaled, csc_matrix(fd.reshape(-1, 1))],
                           [csc_matrix((weights * tangent[:-1]).reshape(1, -1)), csc_matrix([[tangent[-1]]])]], format="csc")
            correction = spsolve(border, -augmented)
            qk += correction[:-1]
            dk += float(correction[-1])
        if not success:
            ds *= 0.5
            if ds < 1e-5:
                rows.append({"status": "corrector_failed", "delta": dk, "ds": ds})
                break
            continue
        old_u, old_delta = u, delta
        u, delta = qk * scales, dk
        row = branch_row(u, delta, intervals, "ok", _iteration + 1, spectral=False)
        row["pseudo_ds"] = step_ds
        row["extended_corrector_residual_inf"] = corrector_residual_inf
        row["extended_corrector_iterations"] = _iteration + 1
        rows.append(row)
        if len(rows) % 50 == 0 or abs(delta - 4.142) < 0.003:
            print("PSEUDO_ARC " + " ".join(f"{k}={v}" for k, v in row.items()))
        new_jac = trapezoid_jacobian(u, delta, intervals)
        fd = trapezoid_parameter_derivative(u, intervals)
        tangent_border = bmat([
            [new_jac @ csc_matrix(np.diag(scales)), csc_matrix(fd.reshape(-1, 1))],
            [csc_matrix((weights * tangent[:-1]).reshape(1, -1)), csc_matrix([[tangent[-1]]])],
        ], format="csc")
        tangent_rhs = np.zeros(u.size + 1)
        tangent_rhs[-1] = 1.0
        candidate = spsolve(tangent_border, tangent_rhs)
        candidate /= np.sqrt(np.dot(weights, candidate[:-1] ** 2) + candidate[-1] ** 2)
        old_tangent_delta = float(tangent[-1])
        if np.dot(weights * candidate[:-1], tangent[:-1]) + candidate[-1] * tangent[-1] < 0.0:
            candidate = -candidate
        if old_tangent_delta * float(candidate[-1]) <= 0.0 and step_ds > 1e-5:
            # Retry the crossing interval with a shorter arc step. This
            # brackets d(delta)/ds=0 without shrinking the whole pre-fold arc.
            rows.pop()
            u, delta = old_u, old_delta
            ds = max(step_ds * 0.5, 1e-5)
            continue
        tangent = candidate
        state_tangent_norm = float(np.sqrt(np.dot(weights, candidate[:-1] ** 2)))
        rows[-1]["tangent_delta"] = float(candidate[-1])
        rows[-1]["state_tangent_norm_scaled"] = state_tangent_norm
        rows[-1]["d_state_d_delta_scaled"] = (
            state_tangent_norm / abs(float(candidate[-1]))
            if abs(float(candidate[-1])) > 1e-14 else float("inf")
        )
        if old_tangent_delta * float(candidate[-1]) <= 0.0:
            # Interpolate the zero of the continuation tangent's parameter
            # component inside the final, locally refined arc-length bracket.
            fraction = old_tangent_delta / (old_tangent_delta - float(candidate[-1]))
            rows[-1]["tangent_zero_fraction_in_last_step"] = fraction
            rows[-1]["fold_delta_tangent_zero"] = old_delta + fraction * (delta - old_delta)
            rows[-1]["fold_delta_bracket_width"] = abs(delta - old_delta)
            rows[-2].update(fixed_scaled_spectrum(old_u, old_delta, intervals))
            rows[-2]["spectrum_point"] = "pseudo_arclength_peak"
            rows[-1].update(fixed_scaled_spectrum(u, delta, intervals))
            rows[-1]["spectrum_point"] = "post_fold"
            rows[-1]["fold_previous_delta"] = old_delta
            rows[-1]["fold_arc_bracket_width"] = step_ds
            rows[-1]["fold_previous_temperature_max_K"] = float(
                300.0 + 15.0 * np.max(old_u.reshape(intervals + 1, 2)[:, 0])
            )
            break
    return rows


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    continuous = collocation_continue()
    hot_mesh = np.linspace(0.0, XI_MAX, 302)
    hot_initial = hot_guess(0.1, hot_mesh).T.reshape(-1)
    hot_u, hot_ok, hot_iterations, hot_residual = newton_correct(hot_initial, 0.1, 301)
    if not hot_ok:
        raise RuntimeError("independent 301-interval hot trapezoid root did not converge")
    hot_profile = hot_u.reshape(302, 2)
    write_csv(HERE / "mf1963_mf_high_hot_trapezoid_301.csv", [
        {"xi": float(xi), "theta": float(theta), "p": float(slope),
         "temperature_K": float(300.0 + 15.0 * theta),
         "Y_A": float(1.0 - GAMMA * theta)}
        for xi, (theta, slope) in zip(hot_mesh, hot_profile)
    ])
    grids = {}
    targets = [0.1, 1.0, 3.0, 4.0, 4.1, 4.12, 4.14, 4.142]
    targets.extend(round(4.142 + i * 0.0001, 7) for i in range(1, 11))
    targets.extend([4.15, 4.2, 4.3, 5.0])
    for intervals in (301, 601, 1201):
        grids[str(intervals)] = discrete_continue(intervals, targets)
    hot_grids = {
        str(intervals): discrete_hot_root(intervals, [0.1, 1.0, 3.0, 4.0, 4.1, 4.12, 4.14, 4.142, 4.143])
        for intervals in (301, 601, 1201)
    }
    arcs = {str(intervals): pseudo_arclength(intervals, points=1200, ds=0.03)
            for intervals in (301, 601, 1201)}
    arc = arcs["301"]
    spectral_by_grid = []
    for intervals in (301, 601, 1201):
        mesh = np.linspace(0.0, XI_MAX, intervals + 1)
        u4, ok4, _, _ = newton_correct(cold_guess(4.0, mesh).T.reshape(-1), 4.0, intervals)
        if ok4:
            spectral_by_grid.append({"point": "delta_4", "intervals": intervals,
                                     **fixed_scaled_spectrum(u4, 4.0, intervals)})
        peak = next((row for row in arcs[str(intervals)]
                     if row.get("spectrum_point") == "pseudo_arclength_peak"), None)
        if peak is not None:
            spectral_by_grid.append({"point": "pseudo_arclength_peak", "intervals": intervals,
                                     **{key: peak[key] for key in (
                                         "delta", "sigma_min_raw", "sigma_max_raw", "condition_raw",
                                         "sigma_min_fixed_scaled", "sigma_max_fixed_scaled",
                                         "condition_fixed_scaled")}})
    arc_mesh = np.linspace(0.0, XI_MAX, 302)
    delta4_seed, delta4_ok, _, _ = newton_correct(
        cold_guess(4.0, arc_mesh).T.reshape(-1), 4.0, 301
    )
    spectral = ([{"spectrum_point": "delta_4_cold_root",
                  **fixed_scaled_spectrum(delta4_seed, 4.0, 301)}]
                if delta4_ok else [])
    spectral.extend(row for row in arc if "spectrum_point" in row)
    report = {
        "model": "theta'=p; p'=p-S(theta,delta); theta(0)=0; p(XI_MAX)=0",
        "xi_max": XI_MAX,
        "beta": BETA,
        "gamma": GAMMA,
        "continuous_collocation": continuous,
        "trapezoid_hot_root_301": {
            "status": "ok", "delta": 0.1, "intervals": 301,
            "iterations": hot_iterations, "residual_inf": hot_residual,
            "temperature_max_K": float(300.0 + 15.0 * np.max(hot_profile[:, 0])),
            "p_inlet": float(hot_profile[0, 1]),
        },
        "trapezoid_by_intervals": grids,
        "trapezoid_hot_branch_by_intervals": hot_grids,
        "trapezoid_pseudo_arclength_301": arc,
        "trapezoid_pseudo_arclength_by_intervals": arcs,
        "spectral_by_grid_fixed_scales": spectral_by_grid,
        "pseudo_arclength_peak_spectrum_fixed_scales": spectral,
        "spectral_scaling": "Fixed scaled Jacobian is diag(row_scale)*J*diag(column_scale): endpoint boundary rows use row_scale=1, trapezoid ODE rows use 1/h, state columns use theta_scale=p_scale=0.05. Raw spectrum uses unscaled J. Pseudo-arclength metric uses the same state scales and delta_scale=1.",
        "note": "Raw sigma_min/condition estimates are discretization- and scaling-dependent; they are diagnostics, not fold proof by themselves.",
    }
    (HERE / "mf1963_mf_high_cold_branch_diagnostic.json").write_text(
        json.dumps(report, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    write_csv(HERE / "mf1963_mf_high_cold_branch_collocation.csv", continuous)
    write_csv(HERE / "mf1963_mf_high_cold_branch_pseudo_arclength.csv", arc)
    write_csv(HERE / "mf1963_mf_high_trapezoid_hot_branch.csv",
              [row for rows in hot_grids.values() for row in rows])
    write_csv(HERE / "mf1963_mf_high_cold_branch_pseudo_arclength_all_grids.csv",
              [row for rows in arcs.values() for row in rows])
    print(f"Wrote independent diagnostic artifacts under {HERE}")


if __name__ == "__main__":
    main()
