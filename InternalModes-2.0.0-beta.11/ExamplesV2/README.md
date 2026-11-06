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
- `waveModeCountsByWavenumber.m` — plot converged, grid-supported and combined leading wave-mode counts at every distinct positive Fourier radius of a specified periodic domain. Uses a fixed physical WKB–Chebyshev grid, two explicit eigenproblem resolutions, and the physical G Gram metric; reports the zero-wavenumber inertial inventory separately.

Run the wave-count example after adding InternalModes and its dependencies to the MATLAB path:

```matlab
addpath ExamplesV2
result = waveModeCountsByWavenumber();
% Or select a new output directory and a finer physical sampling grid:
result = waveModeCountsByWavenumber(Nz=33,outputDirectory="wave-counts-33");
```

The default is a 1 km square periodic domain with 16 × 16 horizontal points, 1 km depth, 25 vertical points, and exponential stratification with surface N = 0.01 s⁻¹ and scale depth 400 m. It requests 24 wave modes (including the external mode) and independently requests 12 inertial modes. Eigenproblem resolutions are 96 and 144; the comparison uses relative equivalent-depth and joint H¹ errors at tolerance 10⁻⁵, and the sampled Gram tolerance is 1%. These are explicit example policies, not universal physical accuracy guarantees. Reaching the candidate ceiling gives a lower bound on usable modes; it does not establish the true maximum. Increasing horizontal resolution extends the examined wavenumber range, whereas changing vertical points changes the fixed physical sampling rule.

Each run writes PNG/PDF figures, wave and inertial count tables, per-mode convergence measurements, per-prefix Gram measurements, the physical quadrature rule and JSON provenance. With no `outputDirectory`, it creates a unique temporary directory and prints its location. Existing output directories are rejected. The source revision and dirty-worktree state are recorded when run from an authoring checkout. Mode convergence and sampled-transform evidence remain separate; neither the plot nor its combined count certifies quadratic products, leakage from arbitrary excluded modes, or absolute solver accuracy. The two-resolution comparison can be repeated at higher resolutions to strengthen the evidence.

See [Resolved-mode convergence and counts by wavenumber](../Documentation/ModeConvergenceAssessment.md) for the measurement definitions, exact count-map contract and interpretation.

The `DiscreteTransform*` scripts isolate individual transform concepts. `DiscreteTransformCustomObjective.m` and `QuadratureWeightRegularizationSweep.m` are investigations of alternative fitting objectives, not recommended production settings.
