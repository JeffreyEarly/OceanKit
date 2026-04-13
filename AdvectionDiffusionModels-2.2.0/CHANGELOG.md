# Version History

## [2.2.0] - 2026-04-13
- simplified the `GriddedStreamfunction` interface and examples, added the new constrained-streamfunction options, and expanded the Site 1 bootstrap workflow
- refactored the gridded-streamfunction bootstrap and trajectory-decomposition antiderivative paths to use the public `BSpline` integration utilities added in `SplineCore 2.2`, and raised the package dependency floor to `SplineCore ^2.2`
- replaced the legacy `kde` and `kde2d` helpers with the shared `KernelDensityEstimate` subsystem, switched bootstrap likelihood and consensus scoring to direct KDE point evaluation, and optimized the one-dimensional and two-dimensional bandwidth-selection paths
- reused observed samples, fast-basis factorizations, and terminal-displacement weights to reduce bootstrap and gridded-projection cost, while simplifying bootstrap summary diagnostics and mesoscale degree-of-freedom bookkeeping
- expanded the estimator website documentation and added three runnable tutorials: stochastic trajectories, the gridded streamfunction fit workflow, and the bootstrap error-assessment and model-choice workflow with rendered model-choice table output

## [2.1.0] - 2026-04-09
- updated the linear-velocity estimation path to the `Distributions` 2.0 named-argument API
- raised the package dependency floor to `Distributions ^2.0`
- migrated the remaining spline call sites to the `SplineCore 2.1` constructor and factory APIs, including the gridded-streamfunction workflows and legacy second-moment helper paths

## [2.0.0] - 2026-03-25
- updated the package for `SplineCore` 2.0 compatibility, including the internal spline call sites and package dependency metadata
- modernized the `Integrators` subsystem with updated constructors, lowerCamel public methods, reorganized class folders, refreshed examples, and expanded inline API documentation
- restored the standard OceanKit website scaffolding with generated class reference pages for the four integrator classes
