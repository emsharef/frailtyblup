"""Scaled normal-map search for the manuscript's constrained BLUP equations.

Scaling is fixed at the projected initial estimate. Every candidate is checked
against the original unscaled equations and normal-cone support functions.
"""

import time
import numpy as np
from scipy.optimize import least_squares
from .derivatives import FastBoundaryProblem


class ScaledBoundaryProblem(FastBoundaryProblem):
    def __init__(self, des, pilot, **kwargs):
        super().__init__(des, **kwargs)
        x = self.project(self.pack(pilot))
        J = self.score_jac(x)
        self.scales = np.ones(self.dim)
        self.scales[: des.p] = 1 / max(1.0, np.max(np.sum(abs(J[: des.p]), axis=1)))
        if not des.subject:
            for k in range(des.p, des.p + 2):
                self.scales[k] = 1 / max(1.0, np.sum(abs(J[k])))
        self.scales[-3:] = 1 / max(1.0, np.max(np.sum(abs(J[-3:]), axis=1)))
        self.original = FastBoundaryProblem(des, **kwargs)

    def normal(self, w):
        x = self.project(w)
        return self.scales * self.score(x) + x - w

    def normal_jac(self, w):
        x = self.project(w)
        P = self.projection_jac(w)
        return (
            self.scales[:, None] * self.score_jac(x) + np.eye(self.dim)
        ) @ P - np.eye(self.dim)

    def certify(self, w):
        x = self.project(w)
        # Translate to the original auxiliary coordinate and apply every original
        # test, including the original unscaled normal-map and support-function tests.
        original_w = x + (w - x) / self.scales
        result = self.original.diagnostics(original_w)
        result["scaled_normal_max"] = float(np.max(abs(self.normal(w))))
        result["translation_shift"] = float(
            np.max(abs(self.original.project(original_w) - x))
        )
        return result

    def solve_scaled(self, pilot, max_nfev, zero=()):
        start = time.monotonic()
        w0 = self.project(self.pack(pilot))
        fixed = self.des.p + np.array(zero, dtype=int)
        w0[fixed] = 0
        free = np.array([k for k in range(self.dim) if k not in fixed])

        def expand(v):
            w = w0.copy()
            w[free] = v
            return w

        opt = least_squares(
            lambda v: self.normal(expand(v))[free],
            w0[free],
            jac=lambda v: self.normal_jac(expand(v))[np.ix_(free, free)],
            xtol=2e-12,
            ftol=None,
            gtol=None,
            max_nfev=max_nfev,
            x_scale="jac",
        )
        w = expand(opt.x)
        x = self.project(w)
        w[fixed] = np.minimum(self.scales[fixed] * self.score(x)[fixed], 0.0)
        result = self.certify(w)
        result.update(
            nfev=opt.nfev,
            njev=opt.njev,
            seconds=time.monotonic() - start,
            solver_message=opt.message,
            solver_success=bool(opt.success),
        )
        return result


def fit(des, initial, max_nfev=300, face_nfev=100, **kwargs):
    pilot = np.r_[initial["hazard"], initial["components_unprojected"]]
    p = ScaledBoundaryProblem(des, pilot, **kwargs)
    projected = p.project(p.pack(pilot))
    attempts = []
    a = p.solve_scaled(pilot, max_nfev)
    a["search"] = "pilot"
    attempts.append(a)
    if not des.subject:
        for zero in [[0, 1], [0], [1]]:
            a = p.solve_scaled(pilot, face_nfev, zero)
            a["search"] = "clinic_zero_" + "".join(map(str, zero))
            attempts.append(a)
    if not any(a["accepted"] for a in attempts):
        a = p.solve_scaled(np.r_[pilot[: des.p], np.zeros(des.e)], max_nfev)
        a["search"] = "zero"
        attempts.append(a)
    for a in attempts:
        a["distance_to_pilot"] = float(
            np.linalg.norm(p.pack(a["estimate"]) - projected)
        )
    accepted = [a for a in attempts if a["accepted"]]
    chosen = min(accepted, key=lambda a: a["distance_to_pilot"]) if accepted else None
    return chosen, attempts, p.scales
