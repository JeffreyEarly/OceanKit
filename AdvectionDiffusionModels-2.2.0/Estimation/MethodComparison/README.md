# Method Comparison

This folder contains apples-to-apples comparison tools for the legacy
Oscroft-style `LinearVelocityField` estimator and the newer
`GriddedStreamfunction` estimator.

The main Site 1 sanity-check path is:

- `buildSite1OscroftFairComparison.m`
- `Examples/FluidsPaperCaseStudy/MakeFigureSite1OscroftFairComparison.m`
- `UnitTests/OscroftFairComparisonUnitTests.m`

## When the two methods line up

For the specific Site 1 fair-comparison settings, the two methods use the same
reduced model family.

Those settings are:

- old method: `EstimateLinearVelocityFieldParameters(..., [u0v0, strain], dof=4)`
- new method: `GriddedStreamfunction.fromTrajectories(..., psiS=[2 2 3], fastS=3, mesoscaleConstraint="zeroVorticity")`
- comparison variables: `uCenter`, `vCenter`, `sigma_n`, `sigma_s`, `sigma`,
  and a branch-aligned principal `theta`

Under these conditions, both methods represent the same centered-frame velocity
field

$$
u(q,r,t) = u_{\mathrm{Center}}(t) + \frac{1}{2}\sigma_n(t) q + \frac{1}{2}\sigma_s(t) r,
$$

$$
v(q,r,t) = v_{\mathrm{Center}}(t) + \frac{1}{2}\sigma_s(t) q - \frac{1}{2}\sigma_n(t) r,
$$

with centered coordinates

$$
q = x - m_x(t), \qquad r = y - m_y(t).
$$

The key point is that the new method does not use a large spatial spline in
this setup. With `psiS=[2 2 3]` and the default `psi` knot construction, there
are no interior `psi` knots, so the mesoscale streamfunction is a single global
tensor-product spline: quadratic in space and cubic in time. Imposing
`mesoscaleConstraint="zeroVorticity"` removes the non-harmonic quadratic terms,
leaving exactly the centered linear zero-vorticity family used by the old
strain-only model.

On the old side, `dof=4` is also a global cubic-in-time basis with four temporal
degrees of freedom, not a four-segment piecewise fit. That is why the overlap
between the two methods is larger than just the constant-in-time case.

So the constant-in-time case is a subset of the overlap, not the whole overlap.
The two methods line up for the full zero-vorticity, globally linear,
cubic-in-time family.

The synthetic comparison test in
`UnitTests/OscroftFairComparisonUnitTests.m` checks exactly this regime and
shows that the old and new summaries agree in common space for `uCenter`,
`vCenter`, `sigma_n`, and `sigma_s`.

## Why the fitted answers can still differ

Even when the model family is the same, the estimators are not the same.

### 1. They construct velocities differently

The old estimator starts from position matrices `x(t), y(t)` and computes
velocities with a second-order finite-difference operator.

The new estimator starts from `TrajectorySpline` objects and uses spline
derivatives to obtain velocities.

On the actual Site 1 trajectories, these two velocity constructions are not
identical. In a local diagnostic on the Site 1 data:

- `uSpline - uFD` had RMS about `1.0e-3 m/s`
- `vSpline - vFD` had RMS about `1.2e-3 m/s`
- pointwise maxima were about `5e-3 m/s`

That difference is large enough to move the fitted strain history.

### 2. They handle the center-of-mass differently

The old estimator uses the sample center of mass directly:

$$
m_x(t) = \mathrm{mean}_k\, x_k(t), \qquad m_y(t) = \mathrm{mean}_k\, y_k(t),
$$

and then finite-differences those time series.

The new estimator first fits a fast spline for the COM trajectory and then
differentiates that spline. On Site 1, the fitted COM position is effectively
identical to the sample mean position, but the COM velocity still differs from
the old finite-difference COM velocity at about `1e-3 m/s` scale.

### 3. They solve different regression problems

The old estimator solves directly for the time-dependent linear coefficients.

The new estimator solves for gauge-reduced streamfunction coefficients, then
recovers `uCenter`, `vCenter`, `sigma_n`, and `sigma_s` from derivatives of the
fitted streamfunction.

So even when both parameterizations span the same set of physical fields, they
do not use the same design matrix, the same state variables, or the same route
from observed trajectories to fitted coefficients.

### 4. The new method is built for asynchronous spline trajectories

The old estimator is based on synchronous position matrices. The new estimator
is built around per-drifter `TrajectorySpline` objects and remains well-defined
for asynchronous sampling.

That broader formulation is part of why the new method is not literally the old
estimator rewritten in different notation.

## What does not explain the difference

The bootstrap scoring does not explain the main residual mismatch in the fair
comparison.

In the staged Site 1 diagnostics, changing from mixed old/new scoring to the
shared common-space score only changes the residual modestly. The main mismatch
is already present after aligning variables and bootstrap resamples, which means
the residual is primarily an estimator difference rather than a score-selection
difference.

The score still matters for which bootstrap member is highlighted in a ranked
envelope, but it does not create the underlying old-versus-new discrepancy.

## Practical interpretation

For the Site 1 fair-comparison settings:

1. the old and new methods can represent the same zero-vorticity linear field
2. they are not numerically the same estimator on real drifter data
3. the main source of mismatch is the upstream treatment of trajectory
   derivatives and COM velocity, not the bootstrap scorer

That is why this folder compares both methods in a shared state space and on
shared bootstrap draws: once those alignments are in place, the remaining
residual can be interpreted as the genuine estimator difference.
