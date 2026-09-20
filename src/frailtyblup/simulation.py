"""Synthetic examples with an exact positive conditional frailty moment law.

This educational generator uses no clinical records or empirical design profiles.
It is separate from the application-informed simulation study in the paper.
"""

import numpy as np


def simulate_demo(
    *,
    n_clusters=30,
    cluster_sizes=(2, 3, 4),
    followup=20.0,
    seed=2026,
    hierarchy="nested",
    alternating=False,
    time_dependent=False,
    curved=True,
    covariance=0.025,
):
    """Return (interval records, truth) for two correlated recurrent event types.

    Cluster effects are independent mean-one two-point variables. Conditional
    subject variances are nu[d] * cluster_effect[d] and conditional cross-type
    covariance is ``covariance``. Bivariate signs construct this law exactly.
    Thinning uses analytic hazard upper bounds; covariate updates occur at fixed
    scheduled times. Alternating episodes have a 0.25-unit washout after type 1.
    Fixed administrative censoring ends follow-up. All times are synthetic.
    """
    if (
        isinstance(n_clusters, bool)
        or not isinstance(n_clusters, int)
        or n_clusters < 1
    ):
        raise ValueError("n_clusters must be a positive integer")
    if hierarchy not in {"nested", "subject"}:
        raise ValueError("invalid hierarchy")
    if not np.isfinite(followup) or followup <= 0:
        raise ValueError("followup must be finite and positive")
    sizes = np.asarray(cluster_sizes)
    if (
        sizes.ndim != 1
        or not len(sizes)
        or np.any(sizes < 1)
        or not np.all(sizes == sizes.astype(int))
    ):
        raise ValueError("cluster_sizes must be positive integers")
    rng = np.random.default_rng(seed)
    sizes = (
        np.ones(n_clusters, dtype=int)
        if hierarchy == "subject"
        else np.resize(sizes.astype(int), n_clusters)
    )
    clinic = np.zeros(2) if hierarchy == "subject" else np.array([0.04, 0.0225])
    nu = np.array([0.10, 0.14])
    beta, gamma = np.array([0.25, -0.20]), np.array([0.20, -0.15])
    base, cap = np.array([0.8, 1.1]), 3.0
    upper = base * (1.32 if curved else 1.0)
    records = []
    subject = 0
    for cluster, size in enumerate(sizes):
        mean = 1 + np.sqrt(clinic) * rng.choice([-1.0, 1.0], 2)
        sd = np.sqrt(nu * mean)
        rho = covariance / np.prod(sd)
        if not np.isfinite(rho) or abs(rho) > 1:
            raise ValueError("covariance is infeasible for this synthetic frailty law")
        for _ in range(size):
            sign0 = rng.choice([-1.0, 1.0])
            sign1 = sign0 if rng.random() < (1 + rho) / 2 else -sign0
            frailty = mean + sd * np.array([sign0, sign1])
            assert np.all(frailty > 0)
            x = float(rng.choice([-1.0, 1.0]))
            z = float(rng.choice([-1.0, 1.0]))
            origin = np.zeros(2)
            time, state, release = 0.0, 0, np.inf
            visit = 4.0 if time_dependent else np.inf
            while time < followup:
                active = (
                    np.ones(2, dtype=bool) if not alternating else np.arange(2) == state
                )
                multiplier = frailty * np.exp(
                    beta * x + (gamma * z if time_dependent else 0)
                )
                rates = multiplier * upper * active
                candidate = (
                    time + rng.exponential() / rates.sum() if rates.sum() else np.inf
                )
                stop = min(candidate, visit, release, followup)
                typ, event = -1, False
                if candidate == stop and stop < followup:
                    typ = int(rng.random() * rates.sum() >= rates[0])
                    age = min(stop - origin[typ], cap)
                    hazard = (
                        base[typ] * (1 + 0.32 * np.sin(1.4 * age + [0.0, 0.7][typ]))
                        if curved
                        else base[typ]
                    )
                    event = rng.random() < hazard / upper[typ]
                for d in (0, 1):
                    if active[d] and stop > time:
                        row = dict(
                            subject=f"s{subject}",
                            cluster=f"c{cluster}",
                            event_type=d,
                            start=time,
                            stop=stop,
                            origin=float(origin[d]),
                            event=int(event and typ == d),
                            x=x,
                        )
                        if time_dependent:
                            row["z"] = z
                        records.append(row)
                if event:
                    if not alternating:
                        origin[typ] = stop
                    elif typ == 0:
                        state, origin[1] = 1, stop
                    else:
                        state, release = 2, stop + 0.25
                if stop == release:
                    state, origin[0], release = 0, stop, np.inf
                if stop == visit:
                    z, visit = -z, visit + 4.0
                time = stop
            subject += 1
    truth = {
        "clinic_variances": clinic.tolist(),
        "residual_covariance": [[nu[0], covariance], [covariance, nu[1]]],
        "coefficients": {
            "x": beta.tolist(),
            **({"z": gamma.tolist()} if time_dependent else {}),
        },
        "age_limits": [cap, cap],
        "seed": seed,
        "n_subjects": subject,
        "n_clusters": n_clusters,
        "hierarchy": hierarchy,
        "baseline": "curved" if curved else "constant",
        "alternating": alternating,
    }
    return records, truth
