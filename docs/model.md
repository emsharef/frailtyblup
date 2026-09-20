# Model and notation

The implementation currently supports two event types, indexed by d = 0, 1. Subjects j belong to independent clusters i. A risk interval carries calendar time t, a specified episode origin, and elapsed age a. The working conditional intensity has the form

$$
Y_{ij}^{(d)}(t)\,U_{ij}^{(d)}\lambda_0^{(d)}(a)\exp\{x_{ij}^{(d)}(t)^T\beta_d\}.
$$

The supplied interval table defines risk Y, age and covariates. The baseline is constant within each of K equally spaced age cells, with separate heights for the two types. The parameterization is unconstrained log-cell heights; the paper's baseline level/shape decomposition describes the same fitted hazard. Cells have the same number K across types but may have different maximum ages.

In the nested model, the cluster frailties have conditional mean one, variances sigma_d squared, and zero cross-type covariance. Given the cluster frailties, the subject frailty means equal the corresponding cluster effects, the conditional subject variances are nu_d squared times the cluster effects, and the within-subject cross-type conditional covariance is theta. Distinct subjects are conditionally independent under the stated model. These are moment assumptions, not a Gaussian or gamma mixing distribution fitted by the software.

The public `clinic_variances` are [sigma_0 squared, sigma_1 squared]. `residual_covariance` has diagonal [nu_0 squared, nu_1 squared] and off-diagonal theta. The internal names `sigma0`, `sigma1`, `nu0`, and `nu1` denote variances, not standard deviations. Do not square them again.

With `hierarchy="subject"`, every subject is an independent unit, clinic variances are zero by construction, and the two subject frailties may still be correlated. The input cluster column is ignored in this mode. Select it only when independence across subjects is appropriate.

For the two types, exposure A is integrated from the supplied covariate segments and step hazards. Counts N and A enter the nested linear prediction equations. The mean equations weight exposure and covariates by those predicted frailties. Normalized second-moment equations update covariance components with the working prediction-error correction. The implementation uses an inverse-free formulation that remains defined when clinic variances vanish or residual covariance loses rank. See Sections 2–3 and Appendices J–L of the [paper](../paper/manuscript.md).

The observed risk history must support the censoring and covariate assumptions of the intended analysis. Interval support for censoring does not establish validity under arbitrary informative dropout. Similarly, coding a time-dependent variable does not establish predictability or the moment conditions required by the theory.
