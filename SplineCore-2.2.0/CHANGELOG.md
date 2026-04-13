# Version History

## [2.2.0] - 2026-04-10
- added public `BSpline.integratedSplineState(...)` and `BSpline.integralMatrixForDataPoints(...)` antiderivative utilities for coefficient-state and interpolation-operator workflows
- refactored `BSpline.cumsum(...)` and tensor-spline integration to reuse the public coefficient-state utility instead of constructing temporary per-slice splines

## [2.1.0] - 2026-04-09
- updated internal distribution construction call sites to the `Distributions` 2.0 named-argument API
- raised the package dependency floor to `Distributions ^2.0`
- refactored `TensorSpline`, `InterpolatingSpline`, and `ConstrainedSpline` so their public constructors are cheap canonical state constructors, while expensive scientific setup moved into explicit factories such as `fromKnotPoints(...)` and `fromGriddedValues(...)`
- added `ConstrainedSpline.fromData(...)` as the preferred one-dimensional fitting factory, while keeping `ConstrainedSpline.fromGriddedValues(...)` as the general rectilinear-grid path
- refactored `TrajectorySpline` to the same constructor/factory persistence model, with `TrajectorySpline(options)` now serving as the canonical solved-state constructor and `TrajectorySpline.fromData(...)` replacing the old positional sample-fitting constructor
- added the public read-only axis vocabulary `knotAxes` and `gridAxes`, backed by the new `SplineAxis` object used for canonical solved-state and persistence paths
- added annotated NetCDF persistence coverage for the spline and constraint classes, including restart paths that preserve solved state without rerunning constrained fitting
- updated examples, README usage, and unit tests to the new constructor and factory model, including dedicated spline persistence regression tests

## [2.0.0] - 2026-03-25
- breaking API cleanup across the core spline classes: public APIs now use spline degree `S` and `knotPoints`, low-level `BSpline` and `TensorSpline` constructors are name-value based, and evaluation is standardized around function-call syntax for values plus `valueAtPoints(..., D=...)` for derivatives
- added tensor-product spline support through `TensorSpline`, including tensor basis construction, tensor-grid evaluation, and tensor transforms such as differentiation, integration, nonlinear power refits, and roots in 1-D
- added noisy-data and constrained fitting through `ConstrainedSpline`, including robust iteratively reweighted least squares, optional correlated-error weighting, local point constraints, global shape constraints, smoothing-matrix analysis, and rectilinear-grid fitting in one or more dimensions
- added dedicated constraint classes `SplineConstraint`, `PointConstraint`, and `GlobalConstraint`, including mask-based point-constraint helpers and monotonicity/positivity constraints
- modernized the implementation with tighter property/state access, more argument validation, simplified constructor behavior, and method organization moved into per-file class methods
- expanded unit-test coverage by splitting tests by class and adding coverage for tensor splines, constrained fits, and constraint objects
- substantially expanded the documentation and website, including generated class reference pages, a class-selection guide, MATLAB comparisons, and tutorial/example coverage for interpolation, basis construction, robust fitting, local constraints, global constraints, and mask-constrained tensor fits

## [1.0.1] - 2026-01-16
- - adding support for out-of-order spline evaluation
