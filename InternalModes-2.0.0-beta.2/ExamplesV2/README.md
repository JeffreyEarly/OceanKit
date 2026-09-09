# V2 examples

These developer-facing scripts exercise the prerelease V2 API. They are not registered as website tutorials.

Start with the scripts that demonstrate normal application workflows:

- `HydrostaticModesExponentialStratification.m` — solve and normalize an aligned hydrostatic family.
- `WaveVortexVerticalGridDesign.m` — choose a point-limited F/G grid and retained mode count, retain negative APV modes, and compare signed projection diagnostics with positive Hilbert-majorant error magnitudes.
- `GeostrophicGeneralizedPotentialEnstrophyModes.m` — compare generalized-potential-enstrophy modes, independently designed mode-root grids, and fitted quadrature weights at several horizontal wavenumbers.
- `GeostrophicModesExponentialStratification.m` — construct generalized-energy APV modes through the public factory.
- `GeostrophicTransformComposition.m` — compose APV and zero-APV coordinates.
- `MeanDensityAnomalyModes.m` — project displacement and synthesize mean pressure.
- `AnalyticalGeostrophicZeroAPVModes.m` — use exact constant/exponential boundary modes.

The `DiscreteTransform*` scripts isolate individual transform concepts. `DiscreteTransformCustomObjective.m` and `QuadratureWeightRegularizationSweep.m` are investigations of alternative fitting objectives, not recommended production settings.
