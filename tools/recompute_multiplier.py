"""Verify all reported normal/t comparisons using the included synthetic records."""

from pathlib import Path
import csv
import gzip
import json
import numpy as np
from scipy.stats import norm, t

ROOT = Path(__file__).resolve().parents[1]


def main():
    old = ROOT / "results/manuscript_v010"
    records = list(map(json.loads, gzip.open(old / "coverage_records.jsonl.gz", "rt")))
    normal_rows = list(csv.DictReader((old / "coverage_summary.csv").open()))
    reported = list(
        csv.DictReader((ROOT / "results/manuscript_v011/reference_coverage.csv").open())
    )
    assert len(normal_rows) == len(reported)
    z = norm.ppf(0.975)
    for row, expected in zip(normal_rows, reported):
        for key, value in row.items():
            assert value == expected[key], (key, value, expected[key])
        selected = [
            r
            for r in records
            if r["case"] == row["case"]
            and r["accepted"]
            and (
                row["stratum"] == "all_accepted"
                or r["boundary"] == (row["stratum"] == "boundary")
            )
        ]
        items = [
            (r, i)
            for r in selected
            for i in r["intervals"]
            if i["method"] == row["method"] and i["parameter"] == row["parameter"]
        ]
        dfs = {r["clusters"] - 1 for r, _ in items}
        assert len(dfs) == 1
        df = dfs.pop()
        critical = t.ppf(0.975, df)
        normal = np.array([abs(i["error"]) <= z * i["se"] for _, i in items])
        student = np.array([abs(i["error"]) <= critical * i["se"] for _, i in items])
        n = len(items)
        assert n == int(row["valid"]) and normal.sum() == int(row["covered"])
        gain = student.astype(int) - normal.astype(int)
        calculated = dict(
            df=df,
            normal_multiplier=z,
            t_multiplier=critical,
            halfwidth_ratio=critical / z,
            t_covered=int(student.sum()),
            t_coverage=student.mean(),
            t_coverage_mcse=np.sqrt(student.mean() * (1 - student.mean()) / n),
            paired_gain=gain.mean(),
            paired_gain_mcse=gain.std(ddof=1) / np.sqrt(n),
        )
        for key, value in calculated.items():
            np.testing.assert_allclose(
                value, float(expected[key]), atol=2e-14, rtol=2e-14
            )
    print(
        json.dumps(
            dict(
                passed=True,
                histories=len(records),
                accepted=sum(r["accepted"] for r in records),
                summary_rows_reproduced=len(reported),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
