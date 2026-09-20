import unittest
import numpy as np
from frailtyblup import from_records


def basic_records():
    return [
        dict(
            subject=f"s{s}",
            cluster="c",
            event_type=d,
            start=0.0,
            stop=1.0,
            origin=0.0,
            event=1,
        )
        for s in range(2)
        for d in (0, 1)
    ] + [
        dict(
            subject=f"s{s}",
            cluster="c",
            event_type=d,
            start=1.0,
            stop=3.0,
            origin=0.0,
            event=0,
        )
        for s in range(2)
        for d in (0, 1)
    ]


class IntervalTests(unittest.TestCase):
    def test_exact_integration_and_left_cell_event(self):
        data = from_records(basic_records(), n_bins=2, age_limits=(2.0, 2.0))
        d = data._design
        hazard = np.log([2.0, 3.0, 5.0, 7.0])
        exposure, derivative = d.exposure(hazard)
        np.testing.assert_allclose(exposure, [[5.0, 12.0], [5.0, 12.0]])
        np.testing.assert_allclose(d.eventsum[:, :, 0], 1.0)
        np.testing.assert_allclose(d.eventsum[:, :, 1], 0.0)
        self.assertEqual(data.metadata["intervals_crossing_age_limit"], 4)

    def test_time_unit_conversion(self):
        data = from_records(
            basic_records(), n_bins=1, age_limits=(3.0, 3.0), time_unit=3.0
        )
        np.testing.assert_allclose(data._design.exposure(np.zeros(2))[0], 1.0)

    def test_overlap_and_missing_fields_rejected(self):
        rows = basic_records()
        with self.assertRaisesRegex(ValueError, "overlapping"):
            from_records(rows + [rows[0]])
        with self.assertRaisesRegex(ValueError, "record 0"):
            from_records([{}])

    def test_subject_membership_and_collinear_covariates(self):
        rows = basic_records()
        rows[-1] = {**rows[-1], "cluster": "other"}
        with self.assertRaisesRegex(ValueError, "different clusters"):
            from_records(rows)
        rows = [{**r, "intercept": 1.0} for r in basic_records()]
        with self.assertRaisesRegex(ValueError, "column rank"):
            from_records(rows, covariates=("intercept",), n_bins=1)

    def test_bad_age_or_type_or_nonfinite(self):
        for key, value in [
            ("origin", 1.0),
            ("event_type", 2),
            ("event", 0.2),
            ("stop", np.nan),
        ]:
            rows = basic_records()
            rows[0] = {**rows[0], key: value}
            with self.subTest(key=key), self.assertRaises(ValueError):
                from_records(rows)

    def test_nested_singletons_not_identifiable(self):
        rows = [{**r, "cluster": r["subject"]} for r in basic_records()]
        with self.assertRaisesRegex(ValueError, "multiple subjects"):
            from_records(rows)

    def test_subject_hierarchy_and_no_event_subjects(self):
        rows = [
            {**r, "event": 0 if r["subject"] == "s1" else r["event"]}
            for r in basic_records()
        ]
        data = from_records(rows, hierarchy="subject", n_bins=1)
        self.assertEqual(data.n_clusters, 2)
        np.testing.assert_equal(data._design.N[1], [0, 0])


if __name__ == "__main__":
    unittest.main()
