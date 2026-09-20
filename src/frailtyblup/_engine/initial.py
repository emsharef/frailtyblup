"""Fit the manuscript's initial estimator and its distinct cluster sandwich."""

import numpy as np
from scipy.optimize import minimize, linprog


def calculate(des, slope=6.0):
    def equations(haz):
        A, Av = des.exposure(haz)
        u = (des.N + 1) / (1 + A)
        rec = (des.eventsum - u[:, :, None] * Av).reshape(des.n, des.p)
        F = np.zeros((des.m, des.p), dtype=rec.dtype)
        np.add.at(F, des.ids, rec)
        return F, u

    def jac(haz):
        J = np.zeros((des.p, des.p))
        for k in range(des.p):
            x = haz.astype(complex)
            x[k] += 1e-24j
            J[:, k] = equations(x)[0].sum(axis=0).imag / 1e-24
        return J

    def value(haz):
        A, _ = des.exposure(haz)
        return float(np.sum((des.N + 1) * np.log1p(A)) - des.eventtotal @ haz)

    C = np.zeros((2 * (des.K - 1), des.p))
    lim = []
    for d in [0, 1]:
        for k in range(des.K - 1):
            ix = d * (des.K - 1) + k
            C[ix, d * des.r + k] = -1
            C[ix, d * des.r + k + 1] = 1
            lim.append(slope * des.caps[d] / des.K)
    lim = np.array(lim)
    haz, pilot = des.pilot()
    constrained = bool(np.any(abs(C @ haz) > lim))
    if constrained:
        opt = minimize(
            lambda h: (value(h), -equations(h)[0].sum(axis=0)),
            haz,
            method="SLSQP",
            jac=True,
            bounds=list(zip(des.lower[: des.p], des.upper[: des.p])),
            constraints=[
                {
                    "type": "ineq",
                    "fun": lambda h: np.r_[lim - C @ h, lim + C @ h],
                    "jac": lambda h: np.r_[-C, C],
                }
            ],
            options={"ftol": 1e-13, "maxiter": 1000},
        )
        haz = opt.x
    for _ in range(8):
        F, u = equations(haz)
        g = F.sum(axis=0)
        if max(abs(g)) < 1e-8:
            break
        step = np.linalg.solve(jac(haz), -g)
        trial = haz + step
        if (
            np.any(abs(C @ trial) > lim)
            or np.any(trial < des.lower[: des.p])
            or np.any(trial > des.upper[: des.p])
        ):
            break
        if value(trial) > value(haz) + 1e-8:
            break
        haz = trial
    F, u = equations(haz)
    g = F.sum(axis=0)
    lp = linprog(
        -g,
        A_ub=np.r_[C, -C],
        b_ub=np.r_[lim, lim],
        bounds=list(zip(des.lower[: des.p], des.upper[: des.p])),
        method="highs",
    )

    if not lp.success:
        raise RuntimeError("Initial linear certification failed: " + lp.message)
    gap = max(0.0, float(g @ (lp.x - haz)))
    active = bool(
        np.min(lim - abs(C @ haz), initial=np.inf) < 1e-6
        or np.min(np.r_[haz - des.lower[: des.p], des.upper[: des.p] - haz]) < 1e-6
    )
    J = jac(haz)
    condition = float(np.linalg.cond(J))
    cov = None
    se = None
    if gap < 1e-4 and not active and condition < 1e12:
        IF = -np.linalg.solve(J, F.T).T
        cov = IF.T @ IF
        se = np.sqrt(np.maximum(0, np.diag(cov)))

    def component_rows(h):
        r = equations(h)[1] - 1
        if des.subject:
            return np.c_[r[:, 0] ** 2, r[:, 1] ** 2, r[:, 0] * r[:, 1]]
        sums = np.zeros((des.m, 2), dtype=r.dtype)
        sq = np.zeros_like(sums)
        cross = np.zeros(des.m, dtype=r.dtype)
        np.add.at(sums, des.ids, r)
        np.add.at(sq, des.ids, r * r)
        np.add.at(cross, des.ids, r[:, 0] * r[:, 1])
        sz = des.sizes[:, None]
        eligible = des.sizes > 1
        if not np.any(eligible):
            raise ValueError(
                "Nested covariance needs at least one nonsingleton cluster"
            )
        # Estimate between-subject covariance from eligible clusters only. Repeat
        # its mean for singleton rows so averaging retains that same estimate.
        B = np.zeros_like(sums)
        B[eligible] = (sums[eligible] ** 2 - sq[eligible]) / (
            sz[eligible] * (sz[eligible] - 1)
        )
        B[~eligible] = B[eligible].mean(axis=0)
        V = sq / sz
        return np.c_[B, V - B, cross / des.sizes]

    comp = component_rows(haz)
    eta = comp.mean(axis=0)
    eta_se = None
    if se is not None and (des.subject or np.all(des.sizes > 1)):
        Dh = np.zeros((des.e, des.p))
        for k in range(des.p):
            h = haz.astype(complex)
            h[k] += 1e-24j
            Dh[:, k] = component_rows(h).mean(axis=0).imag / 1e-24
        EIF = (comp - comp.mean(axis=0)) / des.m + IF @ Dh.T
        eta_se = np.sqrt(np.sum(EIF * EIF, axis=0))
    return {
        "hazard": haz.tolist(),
        "se": None if se is None else se.tolist(),
        "covariance": None if cov is None else cov.tolist(),
        "components_unprojected": eta.tolist(),
        "component_se": None if eta_se is None else eta_se.tolist(),
        "objective_gap_bound": gap,
        "gradient_sum_max": float(max(abs(g))),
        "active_constraint": active,
        "condition": condition,
        "success": gap < 1e-4,
        "labels": des.labels[: des.p],
        "r": des.r,
    }
