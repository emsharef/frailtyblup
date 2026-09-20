"""Validate interval records and integrate step-baseline exposures exactly."""

from dataclasses import dataclass
import csv
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from ._engine.design import ApplicationDesign, add_interval


@dataclass(frozen=True)
class RecurrentData:
    """Prepared bivariate recurrent-event data. Construct with :func:`from_records`.

    ``subject_ids`` and ``cluster_ids`` give the ordering of predictions and
    independent clusters. ``metadata`` reports exposure-window exclusions.
    Arrays in the internal design must not be modified after construction.
    """

    _design: ApplicationDesign
    subject_ids: tuple
    cluster_ids: tuple
    covariates: tuple[str, ...]
    metadata: dict

    @property
    def n_subjects(self) -> int:
        return self._design.n

    @property
    def n_clusters(self) -> int:
        return self._design.m


def from_records(
    records: Iterable[Mapping],
    *,
    covariates: Sequence[str] = (),
    n_bins: int = 3,
    age_limits: Sequence[float] = (3.0, 3.0),
    time_unit: float = 1.0,
    hierarchy: str = "nested",
) -> RecurrentData:
    """Prepare risk intervals for two event types, coded 0 and 1.

    Each record has ``subject``, ``cluster`` (nested hierarchy), ``event_type``,
    ``start``, ``stop``, ``origin`` and ``event``, plus named numeric covariates.
    Risk is observed on (start, stop], event is 0 or 1 at stop, and event age is
    stop - origin. Covariates are constant within each supplied interval.
    Origins must not exceed starts; off-risk periods are omitted, not entered as
    exposure. A terminal event exactly at a cell edge uses the left cell.

    ``age_limits`` and times use the same units. Exposure is truncated at the
    type-specific age limit; endpoint events beyond it are excluded and counted
    in metadata. The hazard's time unit is ``time_unit``. Both types use the same
    number of equal-width bins. Use ``hierarchy='subject'`` for independent
    subjects with correlated type-specific frailties but no cluster effect.
    """
    if (
        isinstance(n_bins, bool)
        or not isinstance(n_bins, (int, np.integer))
        or n_bins < 1
    ):
        raise ValueError("n_bins must be a positive integer")
    limits = np.asarray(age_limits, dtype=float)
    if limits.shape != (2,) or not np.isfinite(limits).all() or np.any(limits <= 0):
        raise ValueError("age_limits must contain two finite positive numbers")
    if not np.isfinite(time_unit) or time_unit <= 0:
        raise ValueError("time_unit must be finite and positive")
    if hierarchy not in {"nested", "subject"}:
        raise ValueError("hierarchy must be 'nested' or 'subject'")
    names = tuple(covariates)
    reserved = {"subject", "cluster", "event_type", "start", "stop", "origin", "event"}
    if len(set(names)) != len(names) or any(
        not isinstance(n, str) or n in reserved for n in names
    ):
        raise ValueError(
            "covariate names must be unique strings distinct from interval columns"
        )
    rows = list(records)
    if not rows:
        raise ValueError("at least one interval is required")
    subject_map, cluster_map, memberships = {}, {}, {}
    segments, spans = [[], []], {}
    retained, excluded = [0, 0], [0, 0]
    truncated = 0
    edges = [np.linspace(0, cap, n_bins + 1) for cap in limits]
    for index, row in enumerate(rows):
        try:
            subject = row["subject"]
            cluster = row["cluster"] if hierarchy == "nested" else subject
            if (
                subject is None
                or cluster is None
                or str(subject) == ""
                or str(cluster) == ""
            ):
                raise ValueError("subject and cluster identifiers must be nonempty")
            for identifier in [subject, cluster]:
                hash(identifier)
                if isinstance(identifier, (float, np.floating)) and not np.isfinite(
                    identifier
                ):
                    raise ValueError("identifiers must be finite")
            typ, event = float(row["event_type"]), float(row["event"])
            if typ not in (0.0, 1.0) or event not in (0.0, 1.0):
                raise ValueError("event_type and event must each be 0 or 1")
            typ, event = int(typ), int(event)
            start, stop, origin = (float(row[k]) for k in ("start", "stop", "origin"))
            cov = np.array([float(row[n]) for n in names])
            if not np.isfinite([start, stop, origin, *cov]).all():
                raise ValueError("times and covariates must be finite")
            if stop <= start or origin > start:
                raise ValueError("require origin <= start < stop")
            if subject in memberships and memberships[subject] != cluster:
                raise ValueError("a subject cannot belong to different clusters")
            memberships[subject] = cluster
            j = subject_map.setdefault(subject, len(subject_map))
            cluster_map.setdefault(cluster, len(cluster_map))
            spans.setdefault((subject, typ), []).append((start, stop))
            retained[typ] += add_interval(
                segments[typ], j, start, stop, origin, cov, event, edges[typ]
            )
            excluded[typ] += int(event and stop - origin > limits[typ])
            truncated += int(stop - origin > limits[typ])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"record {index}: {exc}") from exc
    for key, intervals in spans.items():
        ordered = sorted(intervals)
        if any(b[0] < a[1] for a, b in zip(ordered, ordered[1:])):
            raise ValueError(f"overlapping intervals for subject/type {key!r}")
    ids = np.array([cluster_map[memberships[s]] for s in subject_map])
    sizes = np.bincount(ids)
    if hierarchy == "nested" and not np.any(sizes > 1):
        raise ValueError(
            "nested covariance requires at least one cluster with multiple subjects"
        )
    counts = np.zeros((len(ids), 2))
    for d in (0, 1):
        for j, _, _, _, event in segments[d]:
            counts[j, d] += event
        if retained[d] == 0 or not any(s[2] > 0 for s in segments[d]):
            raise ValueError(
                f"event type {d} requires positive retained exposure and at least one event"
            )
    design = ApplicationDesign(
        ids,
        segments,
        counts,
        list(names),
        n_bins,
        limits,
        time_unit,
        subject=hierarchy == "subject",
        name="user_intervals",
    )
    for d in (0, 1):
        if np.linalg.matrix_rank(design.X[d]) < design.r:
            raise ValueError(
                f"event type {d}: exposure design lacks full column rank; check empty bins, constant or collinear covariates"
            )
    observed = np.zeros(len(ids))
    for d in (0, 1):
        np.add.at(observed, design.rec[d], design.dt[d])
    if np.any(observed == 0):
        raise ValueError("every subject must have retained exposure")
    return RecurrentData(
        design,
        tuple(subject_map),
        tuple(cluster_map),
        names,
        {
            "hierarchy": hierarchy,
            "n_bins": n_bins,
            "age_limits": limits.tolist(),
            "time_unit": float(time_unit),
            "input_intervals": len(rows),
            "retained_events": retained,
            "excluded_events": excluded,
            "intervals_crossing_age_limit": truncated,
        },
    )


def read_csv(path: str | Path, **kwargs) -> RecurrentData:
    """Read the interval schema from CSV, preserving identifier strings."""
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return from_records(csv.DictReader(stream), **kwargs)
