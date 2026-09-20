"""Fit the included, entirely synthetic interval table."""

from pathlib import Path
from frailtyblup import fit, read_csv

data = read_csv(
    Path(__file__).with_name("synthetic_intervals.csv"), covariates=("x",), n_bins=2
)
result = fit(data)
print(result.summary())
if not result.success:
    raise SystemExit("No certified fit; inspect the retained attempt diagnostics.")
