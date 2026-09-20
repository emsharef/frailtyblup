"""Measured covariate changes, unequal clusters including singletons."""

from frailtyblup import fit, from_records
from frailtyblup.simulation import simulate_demo

records, truth = simulate_demo(
    cluster_sizes=(1, 2, 4),
    time_dependent=True,
    followup=25.0,
    n_clusters=40,
    seed=2027,
)
data = from_records(records, covariates=("x", "z"), n_bins=3)
result = fit(data)
print(result.summary())
print("Events excluded by the analysis age window:", data.metadata["excluded_events"])
if not result.success:
    raise SystemExit("No certified fit; inspect the retained attempt diagnostics.")
