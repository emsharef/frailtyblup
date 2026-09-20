# Preparing interval data

One row describes one subject, one event type, and a period of continuous risk with constant covariates. Calendar time and elapsed episode age are different quantities.

| Column | Meaning |
|---|---|
| `subject` | Stable, nonempty subject identifier |
| `cluster` | Independent cluster identifier for the nested model |
| `event_type` | 0 or 1 |
| `start`, `stop` | Calendar-time interval (start, stop], with start < stop |
| `origin` | Calendar-time origin of the current episode; origin <= start |
| `event` | 1 if this type's event occurs at stop, otherwise 0 |
| Named covariate columns | Finite numeric values constant during this interval |

Example: a covariate changes at calendar time 2 during a gap that starts at 0 and ends in an event at 4. Supply two rows with the **same origin 0**: (0, 2], event 0, old value; (2, 4], event 1, new value. Do not reset age when a covariate changes. A change at the event time uses the covariate value just before the event.

For a nonalternating recurrent process, the next episode of that type usually starts at its event time. For alternating onset/end episodes, enter rows only for the type currently at risk, and reset the newly active type's origin at the transition. Omit off-risk washout periods. The package accepts either pattern; it does not guess it from event codes.

At right censoring, include the final risk interval with event 0. Retain subjects with no observed events by providing their exposure intervals. For delayed entry into an existing gap, start exceeds origin. Whether such a left-truncated/prevalent history belongs in the analysis still depends on the observation assumptions in the paper.

## Bins, units, and truncation

`n_bins` specifies equal-width cells between zero and each `age_limits[d]`. Age limits are in the same units as the input times. Exposure beyond a limit is omitted and events beyond it are excluded. `data.metadata` reports these exclusions. This is an explicit analysis-window rule: choose it before looking at fits. Events exactly on a cell boundary are assigned to the left cell.

`time_unit` rescales durations used in fitting. With input days and `time_unit=365.25`, fitted hazard heights are per year; `result.baseline_hazard(ages)` still takes ages in days. The adjacent-log-height slope limit is expressed per chosen time unit. Changing units can therefore change the numerical constraints unless you adjust the slope option accordingly.

Center and scale numerical covariates to useful units. The current coefficient box is [-4, 4]. Do not supply an intercept column: the baseline already supplies the intercept. Constant or collinear columns, empty exposure cells, overlapping intervals within a subject/type, and missing/nonfinite values are rejected. There must be at least one retained event of each type. Imputation, factor encoding and a convention for simultaneous events are analysis decisions outside this API.

The same named covariates appear in both types, with separate fitted coefficients. The internal row representation can express more designs, but the public 0.1 API intentionally has this simpler schema. More than two event types and arbitrary unequal-width bins are not implemented.

## CSV

```python
from frailtyblup import read_csv, fit

data = read_csv("intervals.csv", covariates=("treatment", "age10"),
                n_bins=4, age_limits=(168., 56.), time_unit=28.,
                hierarchy="subject")
result = fit(data)
```

CSV identifiers remain strings. Nested models require at least one cluster with more than one subject to separate cluster and subject components. Singleton clusters are retained: the initial between-subject moment uses eligible clusters, and the final coupled equations include every cluster. This initialization extension is documented in the [algorithm guide](algorithm.md).
