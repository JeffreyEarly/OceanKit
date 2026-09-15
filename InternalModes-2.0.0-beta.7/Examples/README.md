# Internal Modes Examples

This folder contains long-lived package examples and several newer exploratory scripts created while designing the modal-basis and projection API. The newer examples are grouped below by purpose.

Run examples from the package root, for example:

```matlab
cd('/Users/jearly/Documents/OceanKitRepositories/internal-modes')
run('Examples/HydrostaticTransformNumerics.m')
```

## Modal Transform And Projection Examples

These examples exercise the public modal-basis, transform, projection, and spectrum ideas.

| Example | Purpose |
| --- | --- |
| `ProjectionOperatorComparison.m` | Compares canonical weighted-adjoint projections, direct inverses, and weighted pseudoinverses on a hydrostatic mode-adapted grid. This is the compact reproduction of the original quadrature/projection-operator experiment. |
| `FixedGridIGWModeAntialiasing.m` | Mimics the wave-vortex-model situation where hydrostatic modes choose the vertical grid, but nonzero-kappa IGW modes must be projected on that fixed grid. It diagnoses which wave modes are resolvable and which should be filtered. |
| `ObservationalProjectionPartialDepth.m` | Demonstrates an observation-grid projection for partial-depth, irregular samples. It selects a non-contiguous set of resolvable G modes and compares true and recovered coefficients and spectra. |
| `VerticalSpectraAndComponentRoles.m` | Demonstrates the vertical spectrum API and component roles: geostrophic F and G spectra are canonical, wave G is canonical, and wave F is numerical-only. |
| `WaveVortexModelVerticalTransformSketch.m` | Sketches how wave-vortex-model can assemble vertical matrices from generic `InternalModesBasis` and `InternalModesTransform` objects without adding WVM coefficient projectors to `internal-modes`. |
| `CosineHalfWindowEigenmodes.m` | Standalone analytic benchmark for half-window spectra using a finite full-depth cosine basis: first find the visible subspace, then order it with a proxy equivalent-depth operator. |
| `HydrostaticWindowModeComparison.m` | Compares continuous partial-window hydrostatic formulations from the notes: F/HKE, standard G/PE, signed G interior, and an absolute-eta interior variant. It includes a `windowModeNormalization` switch for unit versus equivalent-depth modal scaling, and plots window modes, spectra, full-depth continuations, and spatial-versus-spectral energy checks. |

## Hydrostatic Quadrature Investigations

These scripts are exploratory diagnostics for choosing quadrature nodes and increment weights for hydrostatic F/G transforms. They focus on physical `N2(z_i)`, mode-adapted grids from `G^{N+1}` roots, and Parseval accuracy.

| Example | Purpose |
| --- | --- |
| `EffectiveN2QuadratureComparison.m` | Compares using physical `N2(z_i)` versus an effective `N2` implied by G-mode Christoffel-style weights. This helped separate physical sampling from effective quadrature weights. |
| `PhysicalSamplingQuadratureOptimization.m` | Keeps physical `N2(z_i)` and fixed nodes, then compares geometric increments, F-compatible increments, and joint least-squares increments for F/G Gram accuracy. |
| `PhysicalSamplingChebyshevNodeOptimization.m` | Solves the hydrostatic EVP once, evaluates modes from Chebyshev vectors at trial nodes, and tests whether moving quadrature nodes improves shared physical-sampling weights. |
| `HydrostaticTransformNumerics.m` | Focused hydrostatic-only sanity script for the practical transform workflow. It reports `E_spec`, round-trip errors, condition numbers, band-limited diagnostics, and the constant-stratification DCT/DST limit. |
| `HydrostaticQuadratureWeightInvestigation.m` | Broader fixed-node weight study using `E_spec` as the primary Parseval metric. It compares geometric, F-compatible, normalized LS, regularized LS, band-limited refits, and reweighted LS candidates. |
| `HydrostaticRetainedBandGeneralizedQuadrature.m` | Tests generalized Gaussian-style quadrature for retained modal bands by allowing smooth node motion and refitting positive weights. The main result is that fixed `G^{N+1}` nodes plus retained-band weights capture most of the benefit. |
| `GOnlyBoydWeightsComparison.m` | Tests Boyd/Christoffel-style G-only weights and the analogous F-Christoffel construction. It verifies the constant-stratification DCT/DST limit and compares G-only, F-only, endpoint-split, and coupled F/G weight choices. |

## Practical Reading Order

For the current hydrostatic WVM numerics question, start with:

1. `HydrostaticTransformNumerics.m`
2. `HydrostaticQuadratureWeightInvestigation.m`
3. `GOnlyBoydWeightsComparison.m`

For the broader new transform API, start with:

1. `ProjectionOperatorComparison.m`
2. `FixedGridIGWModeAntialiasing.m`
3. `WaveVortexModelVerticalTransformSketch.m`
4. `VerticalSpectraAndComponentRoles.m`
5. `ObservationalProjectionPartialDepth.m`

## Interpretation Notes

- `E_spec` is the main Parseval-quality diagnostic for a retained modal set. Coefficient round trips can be near machine precision even when the quadrature does not preserve the continuous energy norm.
- The square/native F transform may include an endpoint/Nyquist-like mode needed for interpolation mechanics. Spectrum-safe F diagnostics should use the modes paired with the retained G modes, plus the barotropic mode.
- Constant stratification is a required sanity check: the mode-adapted construction should reduce to DCT/DST-style transforms to numerical precision.
- Node optimization is currently exploratory. The simple practical wins came from choosing better positive increment weights, especially retained-band Gram-fit weights, on fixed mode-adapted nodes.
