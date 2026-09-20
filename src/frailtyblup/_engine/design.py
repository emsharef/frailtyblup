"""Application adapter for the unchanged coupled working-BLUP equations.

Piecewise-constant covariate/age-cell rows; exact exposure integration. Nested
Gaussian *working* elimination is algebra, not a fitted frailty distribution.
"""

import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix


def inv2(a):
    z = np.empty_like(a)
    det = a[..., 0, 0] * a[..., 1, 1] - a[..., 0, 1] * a[..., 1, 0]
    z[..., 0, 0] = a[..., 1, 1] / det
    z[..., 1, 1] = a[..., 0, 0] / det
    z[..., 0, 1] = -a[..., 0, 1] / det
    z[..., 1, 0] = -a[..., 1, 0] / det
    return z


class ApplicationDesign:
    def __init__(
        self, ids, segments, counts, labels, K, caps, unit, subject=False, name=""
    ):
        # segments[d]: (subject row, cell index, exposure length, covariate vector)
        self.ids = np.array(ids, int)
        self.n = len(ids)
        self.m = max(ids) + 1
        self.sizes = np.bincount(ids)
        self.K = K
        self.caps = np.array(caps) / unit
        self.unit = unit
        self.subject = subject
        self.name = name
        self.covlabels = labels
        self.q = len(labels)
        self.r = K + self.q
        self.p = 2 * self.r
        self.e = 3 if subject else 5
        self.N = np.asarray(counts)
        self.X = []
        self.dt = []
        self.rec = []
        self.agg = []
        self.eventsum = np.zeros((self.n, 2, self.r))
        self.labels = [
            f"{x}_{d}"
            for d in [0, 1]
            for x in [*[f"log_cell{k+1}" for k in range(K)], *labels]
        ]
        self.labels += (
            ["nu0", "nu1", "theta"]
            if subject
            else ["sigma0", "sigma1", "nu0", "nu1", "theta"]
        )
        for d in [0, 1]:
            rec = []
            xx = []
            dt = []
            for j, k, length, cov, event in segments[d]:
                v = np.r_[np.eye(K)[k], cov]
                if length > 0:
                    rec.append(j)
                    xx.append(v)
                    dt.append(length / unit)
                if event:
                    self.eventsum[j, d] += event * v
            self.rec.append(np.array(rec))
            self.X.append(np.array(xx))
            self.dt.append(np.array(dt))
            self.agg.append(
                csr_matrix(
                    (np.ones(len(rec)), (rec, np.arange(len(rec)))),
                    shape=(self.n, len(rec)),
                )
            )
        self.eventtotal = self.eventsum.sum(axis=0).reshape(-1)
        self.row_scale = np.r_[np.full(self.p, self.m / self.n), np.ones(self.e)]
        self.lower = np.r_[
            np.tile(np.r_[np.full(K, -12.0), np.full(self.q, -4.0)], 2),
            np.full(self.e - 1, np.log(1e-6)),
            -3.0,
        ]
        self.upper = np.r_[
            np.tile(np.r_[np.full(K, 8.0), np.full(self.q, 4.0)], 2),
            np.full(self.e - 1, np.log(20.0)),
            3.0,
        ]

    def exposure(self, haz):
        dtype = np.result_type(haz, float)
        A = np.zeros((self.n, 2), dtype=dtype)
        Av = np.zeros((self.n, 2, self.r), dtype=dtype)
        for d in [0, 1]:
            w = self.dt[d] * np.exp(self.X[d] @ haz[d * self.r : (d + 1) * self.r])
            A[:, d] = self.agg[d] @ w
            Av[:, d] = self.agg[d] @ (w[:, None] * self.X[d])
        return A, Av

    def pilot(self):
        x = np.zeros(self.p)
        for d in [0, 1]:
            dt = self.dt[d]
            X = self.X[d]
            x[d * self.r : d * self.r + self.K] = np.log(
                (self.N[:, d].sum() + 0.5) / (dt.sum() + 0.5)
            )

        def fn(haz):
            A, Av = self.exposure(haz)
            u = (self.N + 1) / (1 + A)
            val = np.sum((self.N + 1) * np.log1p(A)) - self.eventtotal @ haz
            grad = -(self.eventsum - u[:, :, None] * Av).sum(axis=0).reshape(-1)
            return val / self.n, grad / self.n

        opt = minimize(
            fn,
            x,
            jac=True,
            method="L-BFGS-B",
            bounds=list(zip(self.lower[: self.p], self.upper[: self.p])),
            options={"ftol": 1e-13, "gtol": 1e-9, "maxiter": 1000, "maxls": 40},
        )
        return opt.x, {
            "success": bool(opt.success),
            "gradient_max": float(max(abs(fn(opt.x)[1]))),
            "iterations": int(opt.nit),
            "message": str(opt.message),
        }


def add_interval(out, j, start, stop, origin, cov, event, edges):
    """Split at age-cell boundaries; an endpoint event uses the left cell."""
    lo = max(0.0, start - origin)
    hi = min(stop - origin, edges[-1])
    if hi > lo:
        for k in range(len(edges) - 1):
            length = max(0.0, min(hi, edges[k + 1]) - max(lo, edges[k]))
            if length > 0:
                out.append((j, k, length, cov, 0))
    age = stop - origin
    if event and 0 < age <= edges[-1]:
        k = max(
            0, min(len(edges) - 2, int(np.searchsorted(edges, age, side="left") - 1))
        )
        out.append((j, k, 0.0, cov, event))
        return event
    return 0
