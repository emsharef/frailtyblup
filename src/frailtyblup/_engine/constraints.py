"""Normal-map zero finder for the prespecified constrained thesis score.

Least squares is used only to find a zero of the normal map. A stationary
nonzero residual is a failed search, never a replacement statistical estimate.
"""

import time
import numpy as np
from scipy.optimize import least_squares
from .scores import mean, jacobian, equations, project_residual_covariance


def project_cells(v, limit, lo=-12.0, hi=8.0):
    """Dykstra projection onto the box and adjacent-difference halfspaces."""
    x = np.asarray(v, float).copy()
    K = len(x)
    if np.all((x >= lo) & (x <= hi)) and np.max(abs(np.diff(x)), initial=0) <= limit:
        return x
    # A full box plus one slab for each adjacent pair. Slab projections preserve
    # the pair mean; correction vectors are retained, not alternating clipping.
    corrections = np.zeros((K, K))
    for iteration in range(10000):
        before = x.copy()
        oldcor = corrections.copy()
        y = x + corrections[0]
        x = np.clip(y, lo, hi)
        corrections[0] = y - x
        for k in range(K - 1):
            y = x + corrections[k + 1]
            x = y.copy()
            d = y[k + 1] - y[k]
            change = (d - np.clip(d, -limit, limit)) / 2
            x[k] += change
            x[k + 1] -= change
            corrections[k + 1] = y - x
        if max(np.max(abs(x - before)), np.max(abs(corrections - oldcor))) < 2e-13:
            break
    else:
        raise RuntimeError("Cell projection did not converge")
    return x


class BoundaryProblem:
    def __init__(self, des, clinic_cap=20.0, trace_cap=40.0, slope=6.0):
        self.des = des
        self.clinic_cap = clinic_cap
        self.trace_cap = trace_cap
        self.slope = slope
        self.dim = des.p + des.e
        self.coordinate = np.ones(self.dim)
        self.coordinate[-1] = np.sqrt(2)

    def pack(self, physical):
        return np.asarray(physical) * self.coordinate

    def physical(self, x):
        return np.asarray(x) / self.coordinate

    def project(self, w):
        d = self.des
        x = np.array(w, dtype=float)
        for typ in [0, 1]:
            o = typ * d.r
            x[o : o + d.K] = project_cells(
                x[o : o + d.K], self.slope * d.caps[typ] / d.K
            )
            x[o + d.K : o + d.r] = np.clip(x[o + d.K : o + d.r], -4, 4)
        if not d.subject:
            x[d.p : d.p + 2] = np.clip(x[d.p : d.p + 2], 0, self.clinic_cap)
        x[-3:] = project_residual_covariance(x[-3:], self.trace_cap)
        return x

    def score(self, x):
        return mean(self.des, self.physical(x)) * self.coordinate

    def score_jac(self, x):
        return (
            jacobian(self.des, self.physical(x))
            * self.coordinate[:, None]
            / self.coordinate[None, :]
        )

    def projection_jac(self, w):
        # Exact box and spectral derivatives avoid smoothing across a nearby
        # zero-variance face. The small cell-polytope block uses central differences.
        d = self.des
        P = np.zeros((self.dim, self.dim))
        for typ in [0, 1]:
            o = typ * d.r
            v = w[o : o + d.K]
            limit = self.slope * d.caps[typ] / d.K
            if (
                np.min(v) > -12
                and np.max(v) < 8
                and np.max(abs(np.diff(v)), initial=0) < limit
            ):
                P[o : o + d.K, o : o + d.K] = np.eye(d.K)
            else:
                for k in range(d.K):
                    h = 2e-6 * max(1.0, abs(v[k]))
                    e = np.eye(1, d.K, k).reshape(-1) * h
                    P[o : o + d.K, o + k] = (
                        project_cells(v + e, limit) - project_cells(v - e, limit)
                    ) / (2 * h)
            for k in range(o + d.K, o + d.r):
                P[k, k] = float(-4 < w[k] < 4) if abs(w[k]) != 4 else 0.5
        if not d.subject:
            for k in range(d.p, d.p + 2):
                P[k, k] = (
                    float(0 < w[k] < self.clinic_cap)
                    if w[k] not in (0, self.clinic_cap)
                    else 0.5
                )
        R = np.array([[w[-3], w[-1] / np.sqrt(2)], [w[-1] / np.sqrt(2), w[-2]]])
        lam, V = np.linalg.eigh(R)
        xp = self.project(w)[-3:]
        Q = np.array([[xp[0], xp[2] / np.sqrt(2)], [xp[2] / np.sqrt(2), xp[1]]])
        plam = np.diag(V.T @ Q @ V)
        active = plam > 1e-12
        capped = np.maximum(lam, 0).sum() > self.trace_cap
        diagmap = np.diag(active.astype(float))
        if capped:
            diagmap -= np.outer(active, active) / np.sum(active)
        else:
            for j in range(2):
                if abs(lam[j]) < 1e-14:
                    diagmap[j, j] = 0.5
        if abs(lam[1] - lam[0]) > 1e-12:
            off = (plam[1] - plam[0]) / (lam[1] - lam[0])
        else:
            off = 1.0 if lam[0] > 1e-14 else 0.0 if lam[0] < -1e-14 else 0.5
        for k, E in enumerate(
            [
                np.diag([1.0, 0.0]),
                np.diag([0.0, 1.0]),
                np.array([[0.0, 1.0], [1.0, 0.0]]) / np.sqrt(2),
            ]
        ):
            D = V.T @ E @ V
            L = np.diag(diagmap @ np.diag(D))
            L[0, 1] = L[1, 0] = off * D[0, 1]
            H = V @ L @ V.T
            P[-3:, -3 + k] = [H[0, 0], H[1, 1], np.sqrt(2) * H[0, 1]]
        return P

    def normal(self, w):
        x = self.project(w)
        return self.score(x) + x - w

    def normal_jac(self, w):
        x = self.project(w)
        P = self.projection_jac(w)
        return (self.score_jac(x) + np.eye(self.dim)) @ P - np.eye(self.dim)

    def residual(self, x, tau=1.0):
        return x - self.project(x + tau * self.score(x))

    def diagnostics(self, w):
        d = self.des
        x = self.project(w)
        physical = self.physical(x)
        f = self.score(x)
        physical_score = f / self.coordinate
        R = np.array([[physical[-3], physical[-1]], [physical[-1], physical[-2]]])
        Q = np.array(
            [
                [physical_score[-3], physical_score[-1]],
                [physical_score[-1], physical_score[-2]],
            ]
        )
        trace_active = bool(np.trace(R) > self.trace_cap - 1e-6)
        # With no trace multiplier, these are independent covariance KKT tests.
        eigen = np.linalg.eigvalsh(R)
        qeigen = np.linalg.eigvalsh(Q)
        complement = float(np.trace(R @ Q))
        normal_max = float(np.max(abs(self.normal(w))))
        residuals = {
            str(t): float(np.max(abs(self.residual(x, t)))) for t in [0.1, 1.0, 10.0]
        }
        feasible = float(np.max(abs(self.project(x) - x)))
        eta = physical[d.p :]
        upper = trace_active or (
            not d.subject and bool(np.any(eta[:2] > self.clinic_cap - 1e-6))
        )
        # An independent linear normal-cone test for box/adjacent constraints.
        from scipy.optimize import linprog

        C = []
        limits = []
        for typ in [0, 1]:
            for k in range(d.K - 1):
                row = np.zeros(d.p)
                row[typ * d.r + k] = -1
                row[typ * d.r + k + 1] = 1
                C.extend([row, -row])
                limits.extend([self.slope * d.caps[typ] / d.K] * 2)
        bounds = list(zip(d.lower[: d.p], d.upper[: d.p]))
        opt = linprog(
            -f[: d.p],
            A_ub=np.array(C) if C else None,
            b_ub=np.array(limits) if C else None,
            bounds=bounds,
            method="highs",
        )
        hgap = float(f[: d.p] @ (opt.x - x[: d.p])) if opt.success else float("inf")
        cgap = 0.0
        if not d.subject:
            fg = f[d.p : d.p + 2]
            vg = x[d.p : d.p + 2]
            cgap = float(np.maximum(fg, 0).sum() * self.clinic_cap - fg @ vg)
        # Support function of {R PSD, trace R<=cap} is cap*max(lambda_max(Q),0).
        rgap = float(self.trace_cap * max(0.0, qeigen[-1]) - complement)
        normal_cone_gap = max(0.0, hgap) + max(0.0, cgap) + max(0.0, rgap)
        _, u = equations(d, physical)
        hazard_active = bool(
            np.any(abs(physical[: d.p] - d.lower[: d.p]) < 1e-6)
            or np.any(abs(physical[: d.p] - d.upper[: d.p]) < 1e-6)
        )
        slopes = [
            float(
                np.max(abs(np.diff(physical[typ * d.r : typ * d.r + d.K])), initial=0.0)
                / (d.caps[typ] / d.K)
            )
            for typ in [0, 1]
        ]
        hazard_active = hazard_active or max(slopes) > self.slope - 1e-6
        accepted = (
            normal_max < 2e-8
            and residuals["1.0"] < 2e-8
            and feasible < 1e-9
            and normal_cone_gap < 2e-6
        )
        return dict(
            estimate=physical.tolist(),
            score=physical_score.tolist(),
            normal_score_max=normal_max,
            projected_score_max=residuals,
            feasibility_error=feasible,
            accepted=bool(accepted),
            normal_cone_gap=normal_cone_gap,
            hazard_normal_gap=hgap,
            clinic_normal_gap=cgap,
            residual_normal_gap=rgap,
            residual_eigenvalues=eigen.tolist(),
            residual_score_eigenvalues=qeigen.tolist(),
            residual_complementarity=complement,
            active_upper_covariance=bool(upper),
            active_hazard=bool(hazard_active),
            maximum_slopes=slopes,
            minimum_predictor=float(u.min()),
            hazard_score_max=float(np.max(abs(physical_score[: d.p]))),
        )

    def solve(self, physical_start, max_nfev=300):
        start = time.monotonic()
        x0 = self.project(self.pack(physical_start))
        # Cache repeated projection/score/J evaluations only within each method.
        opt = least_squares(
            self.normal,
            x0,
            jac=self.normal_jac,
            xtol=2e-12,
            ftol=2e-12,
            gtol=2e-12,
            max_nfev=max_nfev,
            x_scale="jac",
        )
        result = self.diagnostics(opt.x)
        result.update(
            auxiliary=opt.x.tolist(),
            start=self.physical(x0).tolist(),
            solver_success=bool(opt.success),
            solver_message=opt.message,
            nfev=opt.nfev,
            njev=opt.njev,
            seconds=time.monotonic() - start,
        )
        return result
