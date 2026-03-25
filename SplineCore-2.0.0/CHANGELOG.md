# Version History

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
