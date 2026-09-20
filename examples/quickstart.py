"""Run after `python -m pip install -e .` from the repository root."""

from frailtyblup import fit, from_records
from frailtyblup.simulation import simulate_demo

records, truth = simulate_demo(seed=2026)
data = from_records(records, covariates=("x",), n_bins=3, age_limits=(3.0, 3.0))
result = fit(data)
print(result.summary())
print("Generating regression coefficients:", truth["coefficients"])
if not result.success:
    raise SystemExit("No certified fit; inspect the retained attempt diagnostics.")
