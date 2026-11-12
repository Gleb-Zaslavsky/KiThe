"""Lift the shared MF-HIGH reduced seed into the six-state numerical system."""

import csv
from pathlib import Path


HERE = Path(__file__).resolve().parent
FIXTURES = HERE.parent / "fixtures"
LAMBDA = 0.1
GAMMA = 0.01
DIFFUSIVITY = 1e-4


def lift(input_name: str, output_name: str) -> None:
    with (FIXTURES / input_name).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pe = float(rows[-1]["xi"])

    with (FIXTURES / output_name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("z", "theta", "q", "C0", "J0", "C1", "J1"))
        for row in rows:
            xi = float(row["xi"])
            theta = float(row["theta"])
            slope = float(row["p"])
            c0 = 1.0 - GAMMA * theta
            j0 = -DIFFUSIVITY * GAMMA * pe * slope
            writer.writerow(
                tuple(
                    f"{value:.17e}"
                    for value in (
                        xi / pe,
                        theta,
                        LAMBDA * pe * slope,
                        c0,
                        j0,
                        1.0 - c0,
                        -j0,
                    )
                )
            )
    print(f"wrote {output_name} ({len(rows)} nodes)")


lift("mf1963_mf_high_two_state_seed.csv", "mf1963_mf_high_six_state_seed.csv")
lift("mf1963_mf_low_two_state_seed.csv", "mf1963_mf_low_six_state_seed.csv")
lift(
    "mf1963_mf_low_uniform_two_state_seed.csv",
    "mf1963_mf_low_uniform_six_state_seed.csv",
)
