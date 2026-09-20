"""Small joint-growth demonstration; all failures are retained in the JSON output."""

import argparse
import json
from pathlib import Path
import numpy as np
from frailtyblup import from_records, fit
from frailtyblup.simulation import simulate_demo


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument(
        "--output", type=Path, default=Path("examples/output/growth.json")
    )
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("replicates must be positive")
    rows = []
    for clusters, followup, bins in [(20, 20.0, 2), (40, 40.0, 3), (80, 80.0, 5)]:
        for rep in range(args.replicates):
            row = dict(
                clusters=clusters,
                followup=followup,
                bins=bins,
                replicate=rep,
                seed=41000 + 100 * clusters + rep,
            )
            try:
                records, truth = simulate_demo(
                    n_clusters=clusters, followup=followup, seed=row["seed"]
                )
                result = fit(from_records(records, covariates=("x",), n_bins=bins))
                row.update(success=result.success, failure=result.failure)
                if result.success:
                    row["regression_error"] = [
                        result.coefficients[d]["x"] - truth["coefficients"]["x"][d]
                        for d in (0, 1)
                    ]
                    row["normal_cone_gap"] = result.selected["normal_cone_gap"]
                row["attempted_searches"] = len(result.attempts)
            except Exception as exc:
                row.update(success=False, failure=type(exc).__name__ + ": " + str(exc))
            rows.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
    for clusters in sorted({r["clusters"] for r in rows}):
        group = [r for r in rows if r["clusters"] == clusters]
        accepted = [r for r in group if r["success"]]
        rmse = (
            np.sqrt(
                np.mean(
                    np.array([r["regression_error"] for r in accepted]) ** 2, axis=0
                )
            )
            if accepted
            else None
        )
        print(
            f"clusters={clusters}: {len(accepted)}/{len(group)} accepted; conditional regression RMSE={rmse}"
        )
    print(
        "Small demonstration only; not a convergence test or reproduction of the paper's study."
    )


if __name__ == "__main__":
    main()
