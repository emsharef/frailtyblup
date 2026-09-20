# Estimation and certification

## Initialization

The initial hazard estimator minimizes the paper's profiled working objective using the predictor (N + 1)/(1 + A), subject to the hazard coefficient box and adjacent-height limits. A linear-program support-function check certifies the initial optimum to the specified numerical tolerance. Its raw component moments are projected to a feasible starting point for the coupled equations.

For nested clusters, the between-subject initial moment averages within-cluster cross-subject products over clusters with size greater than one. Singleton clusters contribute to the subject second moments and the final estimating equations. For the all-nonsingleton designs in the original implementation this is the same initialization. This package extends initialization to mixed singleton/nonsingleton designs and deliberately withholds the old initial component SE calculation for that extension. The final estimating equations are unchanged.

## Coupled search

1. Form the inverse-free nested BLUP and normalized component equations.
2. Fix positive block scalings at the projected initial estimate.
3. Solve the normal map from that start, with analytic score derivatives and a projection derivative.
4. In a nested model, also search the three faces corresponding to one or both zero clinic variances.
5. If no candidate is certified, attempt a zero-covariance start.
6. Select the certified candidate closest to the projected initial estimate in the declared Euclidean/Frobenius coordinates.

Residual off-diagonal covariance is scaled by square root of two in the projection coordinates. The covariance block is constrained to be positive semidefinite with trace at most 40 by default. Clinic variances lie in [0, 20]. Log baseline heights lie in [-12, 8]; regression coefficients lie in [-4, 4]. The default maximum adjacent-log-height slope is 6. Exposed caps and slope are in `FitOptions`; the height and regression boxes are fixed in this release.

SciPy least squares finds a zero of the normal map. It does **not** redefine the statistical estimator as a squared-score minimizer. The solver's success flag alone does not accept a candidate. Original unscaled normal-map and unit-step projected residuals must be below 2e-8, feasibility error below 1e-9, and an independent normal-cone support-function gap below 2e-6. Initial objective certification uses 1e-4. Numerical acceptance can legitimately occur after a solver budget exit if these separate checks pass.

All attempts and solver messages are retained. A global search, uniqueness, and agreement with a particular asymptotic root are not established by this finite search. Computational exceptions still raise exceptions; a completed search without a certified candidate returns `success=False` and no ordinary estimate.

## Predictions and uncertainty

Frailty predictions are linear predictors. They can be negative; they are not clipped. The selected result reports the minimum prediction so this behavior remains visible.

`interior_sandwich()` uses all independent-cluster score rows and the full coupled derivative, including covariance feedback. It raises an error at active hazard/covariance constraints or a near-singular derivative. At an interior fit it returns a diagnostic covariance, whose inferential interpretation still requires the paper's assumptions. It does not assert good small-sample coverage. Boundary confidence intervals, automatic bootstrap inference and a fitted frailty mixing distribution are not part of this release.

## Numerical provenance

The numerical kernels derive from the implementation used for manuscript 0.12. [Kernel provenance](kernel_provenance.json) records original file hashes. Packaging removed filesystem/data dependencies and the unused older interior solver; added validation, configurable search budgets, singleton initialization and one-bin handling; and formatted the code. Tests compare stored original-implementation scores, derivatives and predictions for interior, zero-clinic and singular-residual cases, and independently compare analytic derivatives with complex steps. Tests are numerical regression checks, not proofs of the paper's theorems.
