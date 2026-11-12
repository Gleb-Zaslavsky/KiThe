"""Lightweight contract tests for the independent MF-1963 reference."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path

import numpy as np

import generate_reference as reference
import diagnose_cold_branch as diagnostic


FIXTURES = Path(__file__).resolve().parent / "fixtures"


class ReferenceMathTests(unittest.TestCase):
    def test_source_derivative_matches_centered_difference(self):
        for delta in (1e-4, 0.1):
            for theta in (0.0, 1.0, 50.0, 99.0):
                h = 1e-5
                numerical = (
                    reference.source(theta + h, delta)
                    - reference.source(theta - h, delta)
                ) / (2.0 * h)
                analytic = reference.source_derivative(theta, delta)
                self.assertAlmostEqual(analytic, numerical, delta=2e-8 * max(1.0, abs(analytic)))

    def test_nondimensional_round_trip(self):
        theta = np.array([0.0, 0.5, 20.0, 99.0])
        temperature = reference.T0_K * (1.0 + reference.BETA * theta)
        theta_back = (temperature - reference.T0_K) / reference.T_SCALE
        eta = reference.GAMMA * theta
        reactant = 1.0 - eta
        np.testing.assert_allclose(theta_back, theta, rtol=2e-15, atol=2e-14)
        np.testing.assert_allclose(reactant + eta, 1.0, rtol=0.0, atol=1e-15)
        self.assertAlmostEqual(
            reference.HEAT_RELEASE / reference.CP,
            reference.T0_K * reference.BETA / reference.GAMMA,
        )

    def test_delta_to_arrhenius_round_trip(self):
        for _name, delta, _description in reference.CASES:
            params = reference.dimensional_parameters(delta)
            reconstructed = (
                reference.HEAT_RELEASE
                * reference.LAMBDA
                * reference.RHO0
                * reference.ENERGY_ACTIVATION
                / (reference.CP**2 * reference.MASS_FLUX**2
                   * reference.R_UNIVERSAL * reference.T0_K**2)
                * params["k_at_T0_per_s"]
            )
            self.assertAlmostEqual(reconstructed, delta, delta=2e-15 * delta)
            self.assertGreater(params["k0_per_s"], params["k_at_T0_per_s"])


class ColdBranchContinuationTests(unittest.TestCase):
    def test_trapezoid_parameter_derivative_matches_full_residual_fd(self):
        for intervals in (7, 13):
            nodes = np.arange(intervals + 1, dtype=float)
            theta = 12.0 + 31.0 * (nodes / intervals) ** 1.7 + 4.0 * np.sin(
                2.3 * nodes / intervals
            )
            slope = 0.2 + 1.3 * np.cos(1.9 * nodes / intervals)
            state = np.column_stack((theta, slope)).reshape(-1)
            delta = 1.37
            epsilon = 2e-6

            analytic = diagnostic.trapezoid_parameter_derivative(state, intervals)
            finite_difference = (
                diagnostic.trapezoid_residual(state, delta + epsilon, intervals)
                - diagnostic.trapezoid_residual(state, delta - epsilon, intervals)
            ) / (2.0 * epsilon)
            relative_error = np.linalg.norm(analytic - finite_difference) / max(
                np.linalg.norm(analytic), np.linalg.norm(finite_difference), 1e-14
            )
            max_absolute_error = np.max(np.abs(analytic - finite_difference))

            self.assertEqual(analytic[0], 0.0, "left theta boundary row has no delta term")
            self.assertEqual(analytic[-1], 0.0, "right slope boundary row has no delta term")
            np.testing.assert_array_equal(analytic[1::2], 0.0)
            self.assertGreater(analytic[2], 0.0, "first p-balance derivative must be populated")
            self.assertGreater(analytic[2 * intervals], 0.0,
                               "last p-balance derivative must be populated")
            print(
                f"F_delta FD: intervals={intervals} relative_error={relative_error:.6e} "
                f"max_abs_error={max_absolute_error:.6e}"
            )
            self.assertLess(relative_error, 2e-9)


class FrozenFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        metadata_path = FIXTURES / "metadata.json"
        if not metadata_path.exists():
            raise unittest.SkipTest("run generate_reference.py to create frozen fixtures")
        cls.metadata = json.loads(metadata_path.read_text(encoding="ascii"))

    def test_cases_are_successful_and_physical(self):
        self.assertEqual(len(self.metadata["cases"]), 2)
        for case in self.metadata["cases"]:
            self.assertEqual(case["scipy_status"], 0)
            info = case["diagnostics"]
            self.assertLess(abs(info["bc_theta_inlet"]), 2e-7)
            self.assertLess(abs(info["bc_slope_outlet"]), 2e-7)
            self.assertLess(info["tail_conversion_error"], 2e-7)
            self.assertLess(info["source_integral_identity_abs_error"], 2e-6)
            self.assertGreater(info["source_max"], 0.0)
            self.assertGreater(info["inlet_heat_flux_W_per_m2"], 0.0)
            self.assertTrue(math.isfinite(info["temperature_outlet_K"]))

    def test_profiles_are_finite_monotone_and_bounded(self):
        for case in self.metadata["cases"]:
            with (FIXTURES / case["profile_csv"]).open(newline="", encoding="ascii") as stream:
                rows = list(csv.DictReader(stream))
            self.assertGreater(len(rows), 1000)
            values = np.array([[float(row[key]) for key in (
                "xi", "theta", "dtheta_dxi", "reaction_source", "temperature_K", "Y_A"
            )] for row in rows])
            self.assertTrue(np.isfinite(values).all())
            self.assertTrue(np.all(np.diff(values[:, 0]) > 0.0))
            self.assertGreaterEqual(np.min(np.diff(values[:, 1])), -2e-7)
            self.assertGreaterEqual(values[:, 2].min(), -2e-7)
            self.assertGreaterEqual(values[:, 3].min(), -2e-7)
            self.assertGreaterEqual(values[:, 1].min(), -2e-7)
            self.assertLessEqual(values[:, 1].max(), 100.0 + 2e-7)
            self.assertGreaterEqual(values[:, 5].min(), -2e-7)
            self.assertLessEqual(values[:, 5].max(), 1.0 + 2e-7)
            self.assertAlmostEqual(values[0, 1], 0.0, delta=2e-7)
            self.assertAlmostEqual(values[-1, 2], 0.0, delta=2e-7)

    def test_low_and_high_cases_are_numerically_distinct(self):
        low, high = self.metadata["cases"]
        low_info = low["diagnostics"]
        high_info = high["diagnostics"]
        self.assertGreater(low_info["xi_source_max"], 1000.0)
        self.assertLess(high_info["xi_source_max"], 1.0)
        self.assertGreater(high_info["source_max"], 1000.0 * low_info["source_max"])
        self.assertGreater(high_info["inlet_slope"], 1.0)
        self.assertLess(low_info["inlet_slope"], 1e-3)
        self.assertLess(
            abs(low_info["convective_approx_relative_delta_at_theta_90"]), 0.02
        )

    def test_profile_plots_are_present(self):
        for case in self.metadata["cases"]:
            plot = FIXTURES / case["plot_svg"]
            self.assertTrue(plot.is_file())
            self.assertIn("<svg", plot.read_text(encoding="utf-8")[:1000])


if __name__ == "__main__":
    unittest.main()
