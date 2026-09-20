# FrailtyBLUP

A Python implementation of **Linear Prediction and Frailty Inference for Clustered Recurrent Events**. Fit step baselines, covariate effects, and correlated frailty components using first and second moments rather than a fitted frailty distribution.

[Paper](paper/manuscript.md) · [PDF](paper/paper.pdf) · [Data format](docs/data.md) · [API](docs/api.md) · [Examples](docs/examples.md) · [Reproducibility](docs/reproducibility.md)

## Install

Python 3.10 or later, NumPy and SciPy are required. From this repository:

```bash
python -m pip install .
```

For development, use `python -m pip install -e .`. No package-registry release is required. The working name is `frailtyblup`; availability of that name on package registries has not been established.

## Fit a synthetic example

```python
from frailtyblup import fit, from_records
from frailtyblup.simulation import simulate_demo

records, truth = simulate_demo(seed=2026)
data = from_records(records, covariates=("x",), n_bins=3,
                    age_limits=(3.0, 3.0), hierarchy="nested")
result = fit(data)
print(result.summary())
if result.success:
    print(result.coefficients)
    print(result.residual_covariance)
    result.save("fit.json")
```

The examples use entirely synthetic data. Try `python examples/quickstart.py`, `python examples/time_dependent.py`, `python examples/alternating.py`, and `python examples/from_csv.py`.

## What it handles

- Two recurrent event types with signed within-subject frailty covariance.
- Nested clusters of unequal size, including singleton clusters when other clusters contain multiple subjects, or independent subjects.
- Separate step baselines and regression coefficients for each event type.
- Static and piecewise-constant time-dependent covariates.
- Alternating or nonalternating risk, delayed entry into risk, off-risk periods, and right censoring represented by supplied intervals.
- Zero clinic variances and singular residual covariance through constrained estimating equations.

The software estimates from the risk intervals you supply. It does not infer episode origins, risk states, tie conventions, censoring assumptions, or covariate predictability from an arbitrary table. See the [model](docs/model.md) and [data guide](docs/data.md).

## What a successful fit means

`result.success` means a candidate passed the original, unscaled projected-score, feasibility, and normal-cone checks. All searches remain available in `result.attempts`. Failure is explicit; a failed root search is not relabeled as a fitted estimator. Candidate selection is among the certified candidates found, not a global uniqueness certificate.

The paper proves local long-history results under stated assumptions. Finite-history component bias remains, particularly in sparse clustered designs. Interior sandwich covariance is an explicitly named diagnostic; covariance-boundary confidence intervals are not supplied. [Algorithm and statistical scope](docs/algorithm.md).

## Development

```bash
python -m unittest discover -s tests -v
python -m pip install '.[docs]'
python tools/build_docs.py
```

The generated site is in `site/`. Serve it with `python -m http.server --directory site` to read the paper and documentation. For a source/wheel build, install `build` and run `python -m build`. CI exercises tests, examples, documentation, and an installed wheel.

This is research software, version 0.1.0. [Changes](CHANGELOG.md) · [Contribution guide](CONTRIBUTING.md) · [MIT license](LICENSE).
