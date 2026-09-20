"""Alternating episodes, an off-risk washout, and independent subjects."""

from frailtyblup import fit, from_records
from frailtyblup.simulation import simulate_demo

records, truth = simulate_demo(
    n_clusters=80,
    hierarchy="subject",
    alternating=True,
    time_dependent=True,
    followup=30.0,
    seed=73,
)
data = from_records(records, covariates=("x", "z"), hierarchy="subject", n_bins=3)
result = fit(data)
print(result.summary())
if not result.success:
    raise SystemExit("No certified fit; inspect the retained attempt diagnostics.")
