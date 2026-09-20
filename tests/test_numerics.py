import json
from pathlib import Path
import unittest
import numpy as np
from frailtyblup import FitOptions, fit, from_records
from frailtyblup.simulation import simulate_demo
from frailtyblup._engine.derivatives import evaluate
from frailtyblup._engine.scores import jacobian, equations, project_residual_covariance


class NumericalTests(unittest.TestCase):
    def test_complete_fit_matches_original_manuscript_solver(self):
        reference = json.loads(
            (Path(__file__).parent / "fit_reference.json").read_text()
        )
        rows, _ = simulate_demo(**reference["simulation"])
        result = fit(from_records(rows, covariates=("x",)))
        self.assertTrue(result.success)
        np.testing.assert_allclose(
            result.parameters, reference["estimate"], atol=2e-7, rtol=2e-7
        )
        np.testing.assert_allclose(
            result.initial["hazard"], reference["initial_hazard"], atol=2e-7, rtol=2e-7
        )

    def test_interior_sandwich_includes_all_parameter_blocks(self):
        rows, _ = simulate_demo(
            cluster_sizes=(1, 2, 4),
            time_dependent=True,
            followup=25.0,
            n_clusters=40,
            seed=2027,
        )
        result = fit(from_records(rows, covariates=("x", "z")))
        covariance = result.interior_sandwich()
        self.assertEqual(
            covariance.shape, (len(result.parameters), len(result.parameters))
        )
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-12)
        self.assertGreaterEqual(np.linalg.eigvalsh(covariance).min(), -1e-12)

    def test_frozen_manuscript_scores_and_derivatives(self):
        reference = json.loads((Path(__file__).parent / "reference.json").read_text())
        for case in reference:
            with self.subTest(name=case["name"]):
                records, _ = simulate_demo(**case["simulation"])
                data = from_records(
                    records,
                    covariates=("x",),
                    n_bins=2,
                    hierarchy=case["simulation"]["hierarchy"],
                )
                p = np.array(case["parameters"])
                f, j = evaluate(data._design, p)
                np.testing.assert_allclose(f, case["score"], rtol=2e-12, atol=2e-12)
                np.testing.assert_allclose(j, case["jacobian"], rtol=2e-12, atol=2e-12)
                np.testing.assert_allclose(
                    equations(data._design, p)[1], case["predictions"], atol=2e-12
                )
                np.testing.assert_allclose(
                    j, jacobian(data._design, p), rtol=2e-11, atol=2e-11
                )

    def test_spectral_projection_and_trace_constraint(self):
        for v in [[-1, 2, 0.2], [3, 4, 5], [0, 0, 0]]:
            q = project_residual_covariance(v, trace_cap=2.0)
            r = np.array([[q[0], q[2] / np.sqrt(2)], [q[2] / np.sqrt(2), q[1]]])
            self.assertGreaterEqual(np.linalg.eigvalsh(r).min(), -1e-12)
            self.assertLessEqual(np.trace(r), 2.0 + 1e-12)
            np.testing.assert_allclose(
                project_residual_covariance(q, 2.0), q, atol=1e-12
            )

    def test_one_cell_and_singleton_cluster_fit(self):
        records, _ = simulate_demo(
            n_clusters=20, cluster_sizes=(1, 2, 4), curved=False, followup=30.0, seed=8
        )
        data = from_records(records, covariates=("x",), n_bins=1)
        result = fit(data)
        self.assertTrue(result.success, result.summary())
        self.assertTrue(np.isfinite(result.initial["components_unprojected"]).all())
        self.assertLess(result.selected["normal_cone_gap"], 2e-6)
        self.assertEqual(result.baseline_hazard([0, 1, 3]).shape, (3,))
        json.dumps(result.to_dict(), allow_nan=False)

    def test_success_and_failure_are_distinct(self):
        records, _ = simulate_demo(n_clusters=10, seed=54)
        data = from_records(records, covariates=("x",), n_bins=2)
        result = fit(data, options=FitOptions(max_evaluations=1, face_evaluations=1))
        self.assertFalse(result.success)
        self.assertGreater(len(result.attempts), 0)
        with self.assertRaisesRegex(RuntimeError, "No certified estimate"):
            _ = result.coefficients

    def test_boundary_sandwich_is_not_silently_reported(self):
        records, _ = simulate_demo(seed=2026)
        result = fit(from_records(records, covariates=("x",)))
        self.assertTrue(result.success)
        self.assertLessEqual(result.clinic_variances.min(), 1e-7)
        with self.assertRaisesRegex(ValueError, "boundary"):
            result.interior_sandwich()


if __name__ == "__main__":
    unittest.main()
