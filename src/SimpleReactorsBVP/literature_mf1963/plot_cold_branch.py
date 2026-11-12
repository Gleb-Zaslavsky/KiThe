"""Render MF-02B branch diagnostics from the standalone JSON fixture."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
REPORT = json.loads((HERE / "mf1963_mf_high_cold_branch_diagnostic.json").read_text())


def ok_rows(rows):
    return [row for row in rows if row.get("status") == "ok"]


def plot():
    fig, axes = plt.subplots(3, 1, figsize=(9, 12), constrained_layout=True)
    ax_t, ax_a, ax_s = axes

    for row in REPORT["continuous_collocation"]:
        if row.get("status") == "ok":
            ax_t.plot(row["delta"], row["temperature_max_K"], "o", color="black", ms=3)
            ax_a.plot(row["delta"], row["p_inlet"], "o", color="black", ms=3)

    colors = {"301": "tab:blue", "601": "tab:orange", "1201": "tab:green"}
    for intervals, rows in REPORT["trapezoid_by_intervals"].items():
        good = ok_rows(rows)
        ax_t.plot([r["delta"] for r in good], [r["temperature_max_K"] for r in good],
                  ".-", color=colors[intervals], ms=4, label=f"cold trapezoid {intervals}")
        ax_a.plot([r["delta"] for r in good], [r["p_inlet"] for r in good],
                  ".-", color=colors[intervals], ms=4, label=f"trapezoid {intervals}")

    for intervals, rows in REPORT["trapezoid_hot_branch_by_intervals"].items():
        good = ok_rows(rows)
        ax_t.plot([r["delta"] for r in good], [r["temperature_max_K"] for r in good],
                  "--", color=colors[intervals], alpha=0.8, label=f"hot trapezoid {intervals}")
        ax_a.plot([r["delta"] for r in good], [r["p_inlet"] for r in good],
                  "--", color=colors[intervals], alpha=0.8)

    for intervals, rows in REPORT["trapezoid_pseudo_arclength_by_intervals"].items():
        good = ok_rows(rows)
        ax_t.plot([r["delta"] for r in good], [r["temperature_max_K"] for r in good],
                  color=colors[intervals], lw=1.2, label=f"pseudo-arclength {intervals}")
        ax_a.plot([r["delta"] for r in good], [r["p_inlet"] for r in good],
                  color=colors[intervals], lw=1.2)

    damped = [(0.1, 300.188), (3.0, 308.484), (3.5, 311.364), (4.0, 316.617),
              (4.1, 319.036), (4.142, 321.792)]
    ax_t.plot(*zip(*damped), "s", color="tab:red", label="KiThe Damped cold output")
    ax_t.plot([0.1, 1.0, 4.0, 4.143, 10.0], [1800.0] * 5, "^", color="tab:purple",
              label="hot BVP_sci/CN control")

    ax_t.set_ylabel("maximum temperature (K)")
    ax_t.set_title("MF-HIGH cold and hot branches; fixed domain")
    ax_t.legend(ncol=2, fontsize=8)
    ax_a.set_ylabel("inlet gradient p(0)")
    ax_a.set_title("Cold-branch inlet gradient; dashed curves are hot trapezoid roots")
    ax_a.legend(ncol=2, fontsize=8)

    for row in REPORT["spectral_by_grid_fixed_scales"]:
        if row["point"] == "delta_4":
            ax_s.scatter(row["delta"], row["sigma_min_fixed_scaled"],
                         color=colors[str(row["intervals"])], marker="o")
        else:
            ax_s.scatter(row["delta"], row["sigma_min_fixed_scaled"],
                         color=colors[str(row["intervals"])], marker="x",
                         label=f"fold {row['intervals']} intervals")
    ax_s.set_yscale("log")
    ax_s.set_xlabel("reaction parameter delta")
    ax_s.set_ylabel("sigma_min of fixed-scaled Jacobian")
    ax_s.set_title("Sparse singular-value estimates (circles: delta=4; x: pseudo-arclength peak)")
    ax_s.legend()

    fig.savefig(HERE / "mf1963_mf_high_cold_branch_diagnostics.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    plot()
