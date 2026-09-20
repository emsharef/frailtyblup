"""Analytic derivatives of unchanged inverse-free scores for simulation fitting."""

import numpy as np
from .design import inv2
from .constraints import BoundaryProblem


def evaluate(des, physical):
    n = des.n
    m = des.m
    p = des.p
    e = des.e
    r = des.r
    ids = des.ids
    A = np.zeros((n, 2))
    Av = np.zeros((n, 2, r))
    A2 = np.zeros((n, 2, r, r))
    for d in [0, 1]:
        X = des.X[d]
        w = des.dt[d] * np.exp(X @ physical[d * r : (d + 1) * r])
        A[:, d] = des.agg[d] @ w
        Av[:, d] = des.agg[d] @ (w[:, None] * X)
        A2[:, d] = (
            des.agg[d]
            @ (w[:, None, None] * X[:, :, None] * X[:, None, :]).reshape(len(w), r * r)
        ).reshape(n, r, r)
    eta = physical[p:]
    R = np.array([[eta[-3], eta[-1]], [eta[-1], eta[-2]]])
    I = np.eye(2)
    D = np.zeros((n, 2, 2))
    D[:, 0, 0] = A[:, 0]
    D[:, 1, 1] = A[:, 1]
    B = inv2(I + D @ R)
    T = B @ D
    y = des.N - A
    b = np.einsum("nab,nb->na", B, y)
    if des.subject:
        u = 1 + b @ R.T
        Qr = b[:, :, None] * b[:, None, :] - T
        comp = np.array([Qr[:, 0, 0].mean(), Qr[:, 1, 1].mean(), Qr[:, 0, 1].mean()])
    else:
        G = np.diag(eta[:2])
        t = np.zeros((m, 2, 2))
        s = np.zeros((m, 2))
        np.add.at(t, ids, T)
        np.add.at(s, ids, b)
        C = inv2(I + t @ G)
        a = np.einsum("iab,ib->ia", C, s)
        P = G @ C
        V = C @ t
        mu = a @ G.T
        z = b - np.einsum("nab,nb->na", T, mu[ids])
        u = 1 + mu[ids] + z @ R.T
        Qg = a[:, :, None] * a[:, None, :] - V
        Qr = z[:, :, None] * z[:, None, :] - T + T @ P[ids] @ T
        qr = np.sum(Qr / des.sizes[ids, None, None], axis=0) / m
        comp = np.r_[
            Qg[:, 0, 0].mean(), Qg[:, 1, 1].mean(), qr[0, 0], qr[1, 1], qr[0, 1]
        ]
    score = np.r_[(des.eventsum - u[:, :, None] * Av).sum(axis=0).reshape(-1) / n, comp]
    J = np.zeros((p + e, p + e))
    for k in range(p + e):
        dA = np.zeros((n, 2))
        dR = np.zeros((2, 2))
        dG = np.zeros((2, 2))
        if k < p:
            d = k // r
            kk = k % r
            dA[:, d] = Av[:, d, kk]
        elif k == p + e - 3:
            dR[0, 0] = 1
        elif k == p + e - 2:
            dR[1, 1] = 1
        elif k == p + e - 1:
            dR[0, 1] = dR[1, 0] = 1
        else:
            dG[k - p, k - p] = 1
        dD = np.zeros((n, 2, 2))
        dD[:, 0, 0] = dA[:, 0]
        dD[:, 1, 1] = dA[:, 1]
        dB = -B @ (dD @ R + D @ dR) @ B
        dT = dB @ D + B @ dD
        db = np.einsum("nab,nb->na", dB, y) - np.einsum("nab,nb->na", B, dA)
        if des.subject:
            du = db @ R.T + b @ dR.T
            dQr = db[:, :, None] * b[:, None, :] + b[:, :, None] * db[:, None, :] - dT
            dc = np.array(
                [dQr[:, 0, 0].mean(), dQr[:, 1, 1].mean(), dQr[:, 0, 1].mean()]
            )
        else:
            dt = np.zeros((m, 2, 2))
            ds = np.zeros((m, 2))
            np.add.at(dt, ids, dT)
            np.add.at(ds, ids, db)
            dC = -C @ (dt @ G + t @ dG) @ C
            da = np.einsum("iab,ib->ia", dC, s) + np.einsum("iab,ib->ia", C, ds)
            dP = dG @ C + G @ dC
            dV = dC @ t + C @ dt
            dmu = da @ G.T + a @ dG.T
            dz = (
                db
                - np.einsum("nab,nb->na", dT, mu[ids])
                - np.einsum("nab,nb->na", T, dmu[ids])
            )
            du = dmu[ids] + dz @ R.T + z @ dR.T
            dQg = da[:, :, None] * a[:, None, :] + a[:, :, None] * da[:, None, :] - dV
            dQr = (
                dz[:, :, None] * z[:, None, :]
                + z[:, :, None] * dz[:, None, :]
                - dT
                + dT @ P[ids] @ T
                + T @ dP[ids] @ T
                + T @ P[ids] @ dT
            )
            dr = np.sum(dQr / des.sizes[ids, None, None], axis=0) / m
            dc = np.r_[
                dQg[:, 0, 0].mean(), dQg[:, 1, 1].mean(), dr[0, 0], dr[1, 1], dr[0, 1]
            ]
        dh = -(du[:, :, None] * Av).sum(axis=0)
        if k < p:
            dh[d] -= np.sum(u[:, d, None] * A2[:, d, :, kk], axis=0)
        J[:, k] = np.r_[dh.reshape(-1) / n, dc]
    return score, J


class FastBoundaryProblem(BoundaryProblem):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache = {}

    def _evaluate(self, x):
        if "x" not in self._cache or not np.array_equal(x, self._cache["x"]):
            f, j = evaluate(self.des, self.physical(x))
            self._cache = {
                "x": x.copy(),
                "f": f * self.coordinate,
                "j": j * self.coordinate[:, None] / self.coordinate[None, :],
            }
        return self._cache

    def score(self, x):
        return self._evaluate(x)["f"]

    def score_jac(self, x):
        return self._evaluate(x)["j"]
