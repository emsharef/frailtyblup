"""Inverse-free continuous extension of normalized thesis working-BLUP rows.

No new frailty law, no boundary estimator silently substituted for an old root.
Physical eta order: [sigma0,sigma1,]nu0,nu1,theta (variances, not SDs).
"""

import numpy as np
from .design import inv2


def boundary_moments(A, N, ids, sizes, eta, subject=False):
    n = len(A)
    m = len(sizes)
    dtype = np.result_type(A, eta, float)
    I = np.eye(2, dtype=dtype)
    R = np.array([[eta[-3], eta[-1]], [eta[-1], eta[-2]]], dtype=dtype)
    D = np.zeros((n, 2, 2), dtype=dtype)
    D[:, 0, 0] = A[:, 0]
    D[:, 1, 1] = A[:, 1]
    # B=(I+D R)^-1; T=B D is symmetric and positive semidefinite.
    B = inv2(I + D @ R)
    T = B @ D
    y = N - A
    b = np.einsum("nab,nb->na", B, y)
    if subject:
        Q = b[:, :, None] * b[:, None, :] - T
        return 1 + b @ R.T, np.c_[Q[:, 0, 0], Q[:, 1, 1], Q[:, 0, 1]]
    G = np.diag(eta[:2])
    t = np.zeros((m, 2, 2), dtype=dtype)
    s = np.zeros((m, 2), dtype=dtype)
    np.add.at(t, ids, T)
    np.add.at(s, ids, b)
    C = inv2(I + t @ G)
    a = np.einsum("iab,ib->ia", C, s)
    P = G @ C
    V = C @ t
    mu = a @ G.T
    z = b - np.einsum("nab,nb->na", T, mu[ids])
    Qc = a[:, :, None] * a[:, None, :] - V
    Qr = z[:, :, None] * z[:, None, :] - T + T @ P[ids] @ T
    F = np.zeros((m, 5), dtype=dtype)
    F[:, :2] = np.c_[Qc[:, 0, 0], Qc[:, 1, 1]]
    np.add.at(
        F[:, 2:], ids, np.c_[Qr[:, 0, 0], Qr[:, 1, 1], Qr[:, 0, 1]] / sizes[ids, None]
    )
    return 1 + mu[ids] + z @ R.T, F


def equations(des, physical):
    A, Av = des.exposure(physical[: des.p])
    u, comp = boundary_moments(
        A, des.N, des.ids, des.sizes, physical[des.p :], des.subject
    )
    rec = (des.eventsum - u[:, :, None] * Av).reshape(des.n, des.p)
    F = np.zeros((des.m, des.p + des.e), dtype=rec.dtype)
    np.add.at(F[:, : des.p], des.ids, rec)
    F[:, des.p :] = comp
    return F, u


def mean(des, physical):
    return equations(des, physical)[0].mean(axis=0) * des.row_scale


def jacobian(des, physical):
    J = np.empty((len(physical), len(physical)))
    for k in range(len(physical)):
        v = np.asarray(physical, dtype=complex)
        v[k] += 1e-24j
        J[:, k] = mean(des, v).imag / 1e-24
    return J


def project_residual_covariance(v, trace_cap=40.0):
    """Euclidean projection in (R00,R11,sqrt(2)*R01) coordinates.

    The trace cap is an explicit compactness constraint, not a variance prior.
    """
    R = np.array([[v[0], v[2] / np.sqrt(2)], [v[2] / np.sqrt(2), v[1]]])
    w, V = np.linalg.eigh(R)
    positive = np.maximum(w, 0.0)
    if positive.sum() > trace_cap:
        u = np.sort(w)[::-1]
        css = np.cumsum(u)
        rho = np.flatnonzero(u - (css - trace_cap) / np.arange(1, len(w) + 1) > 0)[-1]
        level = (css[rho] - trace_cap) / (rho + 1)
        positive = np.maximum(w - level, 0.0)
    P = (V * positive) @ V.T
    return np.array([P[0, 0], P[1, 1], np.sqrt(2) * P[0, 1]])
