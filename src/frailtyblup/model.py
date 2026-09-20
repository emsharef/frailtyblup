"""Public fitting interface and explicit numerical certification results."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from .data import RecurrentData
from ._engine.initial import calculate
from ._engine.solver import fit as solve
from ._engine.constraints import BoundaryProblem
from ._engine.scores import equations
from ._engine.derivatives import evaluate


@dataclass(frozen=True)
class FitOptions:
    """Search budgets and compact parameter constraints used in the paper.

    Variance caps are variances, not standard deviations. The slope constraint
    applies to adjacent log-baseline heights divided by cell width, measured in
    the chosen time unit. Log heights are bounded by [-12, 8] and regression
    coefficients by [-4, 4]; center/scale covariates before fitting.
    """

    max_evaluations: int = 300
    face_evaluations: int = 100
    clinic_variance_cap: float = 20.0
    residual_trace_cap: float = 40.0
    max_log_hazard_slope: float = 6.0

    def __post_init__(self):
        for key in ("max_evaluations", "face_evaluations"):
            v = getattr(self, key)
            if isinstance(v, bool) or not isinstance(v, int) or v < 1:
                raise ValueError(f"{key} must be a positive integer")
        for key in (
            "clinic_variance_cap",
            "residual_trace_cap",
            "max_log_hazard_slope",
        ):
            if not np.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} must be finite and positive")


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


@dataclass
class FitResult:
    """A selected, certified candidate or an explicit unsuccessful fit.

    ``success`` means the unscaled projected-score and normal-cone checks passed.
    It does not certify a global root search, uniqueness, or statistical accuracy.
    All attempts are retained; an unsuccessful search is never returned as an
    ordinary estimate. Frailty predictions may be negative: the few-moment linear
    predictor is not clipped or replaced by a posterior expectation.
    """

    data: RecurrentData
    options: FitOptions
    initial: dict
    selected: dict | None
    attempts: list
    failure: str | None = None

    @property
    def success(self) -> bool:
        return self.selected is not None

    def _require(self):
        if not self.success:
            raise RuntimeError(
                f"No certified estimate: {self.failure}; inspect initial and attempts"
            )

    @property
    def parameters(self) -> np.ndarray:
        self._require()
        return np.array(self.selected["estimate"])

    @property
    def coefficients(self) -> dict:
        """Regression coefficients by event type (baseline heights excluded)."""
        d, p = self.data._design, self.parameters
        return {
            typ: dict(
                zip(
                    self.data.covariates,
                    map(float, p[typ * d.r + d.K : (typ + 1) * d.r]),
                )
            )
            for typ in (0, 1)
        }

    @property
    def residual_covariance(self) -> np.ndarray:
        p = self.parameters
        return np.array([[p[-3], p[-1]], [p[-1], p[-2]]])

    @property
    def clinic_variances(self) -> np.ndarray:
        d, p = self.data._design, self.parameters
        return np.zeros(2) if d.subject else p[d.p : d.p + 2]

    @property
    def frailty_predictions(self) -> np.ndarray:
        """Subject-by-type predictions in ``data.subject_ids`` order."""
        return equations(self.data._design, self.parameters)[1]

    def baseline_hazard(self, ages, event_type: int = 0) -> np.ndarray:
        """Step heights per time_unit; ages are in the input's original units.

        Values at interior cell edges use the left cell. No extrapolation beyond
        the analysis age window is provided.
        """
        if event_type not in (0, 1):
            raise ValueError("event_type must be 0 or 1")
        d, p = self.data._design, self.parameters
        ages = np.asarray(ages, dtype=float)
        cap = self.data.metadata["age_limits"][event_type]
        if not np.isfinite(ages).all() or np.any(ages < 0) or np.any(ages > cap):
            raise ValueError("ages must be finite and within the fitted age window")
        bins = np.clip(
            np.searchsorted(np.linspace(0, cap, d.K + 1), ages, side="left") - 1,
            0,
            d.K - 1,
        )
        return np.exp(p[event_type * d.r : event_type * d.r + d.K])[bins]

    def interior_sandwich(self) -> np.ndarray:
        """Full cluster sandwich diagnostic, available only at interior fits.

        Returns physical-coordinate covariance in ``parameter_names`` order.
        Its asymptotic justification requires the paper's long-history and
        consistency assumptions; this method does not certify those assumptions.
        Boundary fits raise ValueError rather than receiving interior intervals.
        """
        self._require()
        d = self.data._design
        if (
            self.selected["active_hazard"]
            or self.selected["active_upper_covariance"]
            or np.linalg.eigvalsh(self.residual_covariance).min() <= 1e-7
            or (not d.subject and self.clinic_variances.min() <= 1e-7)
        ):
            raise ValueError(
                "interior sandwich is unavailable at an active constraint/covariance boundary"
            )
        f = equations(d, self.parameters)[0] * d.row_scale
        j = evaluate(d, self.parameters)[1]
        if d.m <= 1 or np.linalg.cond(j) >= 1e12:
            raise ValueError(
                "insufficient clusters or ill-conditioned score derivative"
            )
        influence = -np.linalg.solve(j, (f - f.mean(axis=0)).T).T
        return influence.T @ influence / d.m**2

    @property
    def parameter_names(self) -> tuple:
        return tuple(self.data._design.labels)

    def to_dict(self) -> dict:
        """JSON-ready estimates and diagnostics; does not export input records."""
        return _jsonable(
            {
                "success": self.success,
                "failure": self.failure,
                "parameter_names": self.parameter_names,
                "initial": self.initial,
                "selected": self.selected,
                "attempts": self.attempts,
                "options": vars(self.options),
                "data": self.data.metadata,
                "n_subjects": self.data.n_subjects,
                "n_clusters": self.data.n_clusters,
            }
        )

    def save(self, path: str | Path):
        """Write estimates and all search diagnostics as standard JSON."""
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, allow_nan=False) + "\n"
        )

    def summary(self) -> str:
        lines = [
            f"frailtyblup: {'certified candidate' if self.success else 'unsuccessful search'}",
            f"Subjects: {self.data.n_subjects}; independent clusters: {self.data.n_clusters}",
        ]
        if self.success:
            lines.extend(
                [
                    f"Regression: {self.coefficients}",
                    f"Clinic variances: {self.clinic_variances}",
                    f"Residual covariance:\n{self.residual_covariance}",
                    f"Normal-cone gap: {self.selected['normal_cone_gap']:.3g}",
                    "Numerical certification is not a finite-sample accuracy or inference guarantee.",
                ]
            )
        else:
            lines.append(
                f"Failure: {self.failure}; {len(self.attempts)} search attempts recorded"
            )
        return "\n".join(lines)


def fit(data: RecurrentData, *, options: FitOptions | None = None) -> FitResult:
    """Fit the initial predictor then the jointly constrained BLUP equations.

    Uses the manuscript's fixed positive score scaling, normal-map searches and
    nearest-initial-candidate selection among certified candidates found. Nested
    models search the three zero-clinic faces as well as the full start. A zero
    covariance start is attempted if no preceding candidate is accepted.
    """
    if not isinstance(data, RecurrentData):
        raise TypeError("data must be prepared with from_records or read_csv")
    options = options or FitOptions()
    if not isinstance(options, FitOptions):
        raise TypeError("options must be FitOptions")
    d = data._design
    initial = calculate(d, slope=options.max_log_hazard_slope)
    kwargs = dict(
        clinic_cap=options.clinic_variance_cap,
        trace_cap=options.residual_trace_cap,
        slope=options.max_log_hazard_slope,
    )
    problem = BoundaryProblem(d, **kwargs)
    pilot = problem.pack(np.r_[initial["hazard"], initial["components_unprojected"]])
    feasible = np.max(abs(problem.project(pilot)[: d.p] - pilot[: d.p])) < 1e-8
    if not initial["success"] or not feasible:
        return FitResult(
            data, options, initial, None, [], "initial estimator failed certification"
        )
    selected, attempts, _ = solve(
        d,
        initial,
        max_nfev=options.max_evaluations,
        face_nfev=options.face_evaluations,
        **kwargs,
    )
    return FitResult(
        data,
        options,
        initial,
        selected,
        attempts,
        None if selected is not None else "no coupled candidate passed certification",
    )
