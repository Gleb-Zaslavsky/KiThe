"""Create shared two-state MF meshes/seeds from the frozen references."""

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
CASES = (("mf_high", 301), ("mf_low", 2049))


def interpolate(points: list[tuple[float, float]], x: float) -> float:
    upper = next((i for i, point in enumerate(points) if point[0] >= x), len(points) - 1)
    upper = max(1, upper)
    x0, y0 = points[upper - 1]
    x1, y1 = points[upper]
    fraction = min(1.0, max(0.0, (x - x0) / (x1 - x0)))
    return y0 + fraction * (y1 - y0)


for case_name, node_count in CASES:
    source = HERE.parent / "fixtures" / f"mf1963_{case_name}.csv"
    output = HERE.parent / "fixtures" / f"mf1963_{case_name}_two_state_seed.csv"
    with source.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    theta = [(float(row["xi"]), float(row["theta"])) for row in rows]
    slope = [(float(row["xi"]), float(row["dtheta_dxi"])) for row in rows]
    xi_max = theta[-1][0]

    if case_name == "mf_low":
        source_rows = [
            (float(row["xi"]), float(row["theta"]), float(row["dtheta_dxi"]))
            for row in rows
        ]
        seed_rows = []
        for index in range(node_count):
            position = index * (len(source_rows) - 1) / (node_count - 1)
            lower = int(position)
            upper = min(lower + 1, len(source_rows) - 1)
            fraction = position - lower
            seed_rows.append(
                tuple(
                    source_rows[lower][column]
                    + fraction * (source_rows[upper][column] - source_rows[lower][column])
                    for column in range(3)
                )
            )
    else:
        seed_rows = [
            (
                xi := xi_max * index / (node_count - 1),
                interpolate(theta, xi),
                interpolate(slope, xi),
            )
            for index in range(node_count)
        ]

    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("xi", "theta", "p"))
        for seed_row in seed_rows:
            writer.writerow(
                tuple(f"{value:.17e}" for value in seed_row)
            )
    print(f"wrote {output} ({len(seed_rows)} nodes, xi_max={xi_max:.17e})")

with (HERE.parent / "fixtures" / "mf1963_mf_low.csv").open(
    newline="", encoding="utf-8"
) as stream:
    low_rows = list(csv.DictReader(stream))
low_theta = [(float(row["xi"]), float(row["theta"])) for row in low_rows]
low_slope = [(float(row["xi"]), float(row["dtheta_dxi"])) for row in low_rows]
low_xi_max = low_theta[-1][0]
uniform_output = HERE.parent / "fixtures" / "mf1963_mf_low_uniform_two_state_seed.csv"
with uniform_output.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(("xi", "theta", "p"))
    for index in range(2049):
        xi = low_xi_max * index / 2048
        writer.writerow(
            tuple(
                f"{value:.17e}"
                for value in (
                    xi,
                    interpolate(low_theta, xi),
                    interpolate(low_slope, xi),
                )
            )
        )
print(f"wrote {uniform_output} (2049 uniform nodes)")
