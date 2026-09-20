# Python API

## Data

`from_records(records, *, covariates=(), n_bins=3, age_limits=(3., 3.), time_unit=1., hierarchy="nested")`

Returns `RecurrentData`. Its public attributes include `subject_ids`, `cluster_ids`, `covariates`, `metadata`, `n_subjects`, and `n_clusters`. Subject ordering follows first appearance in the input; predictions use that order. See the [interval schema](data.md). Internal arrays are implementation details and must not be changed.

`read_csv(path, **kwargs)` uses the same arguments and validation, reading interval records from a UTF-8 CSV file.

## Fit

`fit(data, *, options=None)` returns `FitResult`. `options` is `FitOptions(max_evaluations=300, face_evaluations=100, clinic_variance_cap=20., residual_trace_cap=40., max_log_hazard_slope=6.)`.

| Result attribute/method | Meaning |
|---|---|
| `success`, `failure` | Numerical acceptance or reason for an unsuccessful completed search |
| `initial` | Initial hazard fit, component moments and certification diagnostics |
| `attempts` | All searches, accepted and failed, including estimates and diagnostics |
| `selected` | The chosen certified candidate, or None |
| `coefficients` | Dictionary from event type to named regression coefficients |
| `clinic_variances` | Length-two array of cluster variances; zero for subject hierarchy |
| `residual_covariance` | Symmetric two-by-two subject covariance component matrix |
| `parameters`, `parameter_names` | Full physical parameter vector and its ordering |
| `frailty_predictions` | Subject-by-type linear predictions |
| `baseline_hazard(ages, event_type=0)` | Step heights at ages in original input units, rates per time_unit |
| `interior_sandwich()` | Full physical-coordinate cluster covariance diagnostic, subject to interior checks |
| `summary()` | Human-readable fit summary |
| `to_dict()`, `save(path)` | JSON-ready estimates/diagnostics and file serialization |

Accessing an estimate on an unsuccessful result raises `RuntimeError`. Nonfinite diagnostic scalars are represented as JSON null; input records and individual identifiers are not serialized by `save`. `frailty_predictions` is available separately and is not included automatically.

## Synthetic generator

`frailtyblup.simulation.simulate_demo(...)` returns interval records and a truth dictionary. Useful arguments are `n_clusters`, `cluster_sizes`, `followup`, `seed`, `hierarchy`, `alternating`, `time_dependent`, `curved`, and `covariance`. The bounded positive frailty generator satisfies the specified conditional first two moments. No clinical data or empirical patient profiles are used. It is a teaching/example generator, not the historical simulation engine behind the paper's entire study.
