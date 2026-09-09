# Version History

## Unreleased

## [2.0.0-beta.3] - 2026-09-08

- Scale Chebyshev coefficient unknowns in the V2 spectral eigensolve to reduce derivative-amplified roundoff in long external waves. Recover the original native coefficients before finite-mode filtering and normalization; retain the physical EVP, boundary conditions, requested modes and finite-difference behavior. Add direct and sampled vertical-momentum regressions across coordinate choices and spectral resolutions.

## [2.0.0-beta.2] - 2026-09-08

- Added `IMBasisCollection` and `IMProjection` with hidden basis-owned pairing descriptors, structured assessment results, and separate explicit checking functions. Scientific labels, unequal page counts, supplied samples/weights, signed physical duals, positive error norms, and unsupported capabilities remain explicit.
- Added exact bulk wave construction with shared preparation, ordered/repeated wavenumber mapping, independent zero-wavenumber inertial counts, and bounded continuous evaluation. The shared scalar/bulk solve preserves equation equilibration and consistent finite-mode filtering.
- Added `IMProductAssessmentPlan`, factor/output structures, budget reservation before evaluation, independent-reference diagnostics, coverage comparison, and advisory decisions that preserve scientific counts and endpoints.
- Preserved direct-projection and prefix provenance. The three new public classes replace the unreleased recipe, inventory, and assessment-result classes; released transform APIs remain available.

## V2 development through 2.0.0-beta.1
- stabilized the analytical negative-APV characteristic determinant for strongly localized endpoint modes by cancelling known Bessel exponential factors before root finding; added independent spectral and endpoint-residual checks
- equilibrated generalized eigenproblem rows before solving, reducing high-resolution APV mode and endpoint-response errors while preserving the physical matrices for mode classification and diagnostics
- made spectral stretched-coordinate maps, inverse maps, Jacobians, and second derivatives use one smooth Chebfun representation, removing the inconsistent trapezoidal/PCHIP and endpoint-gradient error floor in WKB and density coordinates
- corrected solved G-mode normalization and projection to use the canonical EVP interior weight, including `(N2-f0^2)/g` for fixed-wavenumber waves and `(N2-omega^2)/g` for fixed-frequency waves; positive majorants now also take the absolute interior weight
- removed the obsolete legacy fixed-frequency normalization path from production solvers, analytical helpers, examples, tests, and user documentation
- added the induced Hilbert-majorant APIs `majorantInnerProduct`, `majorantGramMatrix`, `majorantNorm`, and `targetMajorantGramMatrix`; retained signed Pontryagin projection and spectra; and enabled coupled APV quadratic certification with retained negative modes by measuring relative errors in the positive majorant
- moved singular-pencil rejection into generic eigenpair processing and made zero-norm classification mode-local, preventing near-null generalized eigenvectors and unrelated large columns from corrupting physical mode selection
- identified the integrated V2 development package as `2.0.0-beta.1`; this remains a prerelease and does not declare or publish a stable V2 API
- reconciled the scalar and aligned transform contracts, centralized retained-band fitting and normalization policy, simplified redundant validation, and qualified the V2 API, examples, generated reference, and installed package as one integration increment
- replaced the V2 analytical SQG-profile and availability-report APIs with concrete constant- and exponential-stratification `geostrophicZeroAPVModesAtWavenumber` constructors, exact canonical `F/G` endpoint responses, and an `IMAnalyticalGeostrophicZeroAPVModesBasis` supporting the numerical quadratic-form and rotation contracts
- replaced the projected surface-geostrophic API with coefficient-independent `IMGeostrophicZeroAPVModes` boundary-response problems, page-shaped canonical `F/G` bases, explicit physical-energy and endpoint-buoyancy forms, and boundary-depth, surface-buoyancy, and custom symmetric-pencil rotations
- added the public v2 `IMInternalModes.geostrophicAPVModes` generalized-energy APV factory, signed endpoint-limit mappings and metadata, exact numerical zero retention, volume-only `depth` normalization, and exact constant-stratification trigonometric/affine/hyperbolic branches with independent factory, endpoint, normalization, analytical-root, and refined-spectral regression coverage
- added exact exponential-stratification generalized-energy geostrophic APV modes, including positive, zero, and negative equivalent-depth branches, free-surface and rigid-lid conventions, all finite/zero/infinite endpoint limits, signed-spectrum diagnostics, and analytical-versus-numerical regression coverage
- added v2 scalar discrete Galerkin transforms on caller-supplied points and quadrature weights, including forward and inverse matrices, forward and back transform methods, sampled metrics, continuous Gram targets, and transform-quality diagnostics
- added fixed-point v2 quadrature fitting with normalized Gram objectives, nonnegative full-depth constraints, custom least-squares systems, and fitted-versus-geometric diagnostics
- made fixed-point fitting depth-scale robust and added v2 walkthroughs for default and custom scalar discrete-transform workflows
- added v2 `pointsFromModeRoots` grids with automatic auxiliary-mode solves and native Chebyshev root finding for all spectral coordinate choices
- compressed the normalized Gram Frobenius fit to independent symmetric mode pairs
- added exact point-limited scalar transform construction and `IMDiscreteTransformAssessment`, including one-rule retained-prefix Gram diagnostics, optional rejected-mode leakage and scalar quadratic-aliasing policies, strict explicit-band failures, and constant/exponential/DCT-I regression coverage
- added aligned internal-mode `F/G` discrete transforms with one shared point-and-weight rule, variable-qualified synthesis and projection matrices, zero-column active projectors, endpoint and physical metadata snapshots, stacked shared-weight fitting, worst-channel Gram and leakage policies, and coupled `FF->F`, `FG->G`, and `GG->F` quadratic-aliasing assessment
- added `IMGeostrophicTransform` to compose generalized-energy APV transforms with boundary-normalized or rotated geostrophic zero-APV coordinates, including admissible-state inversion, generic volume-plus-endpoint source projection, explicit wavenumber pages, singularity diagnostics, and arbitrary trailing field dimensions
- added generalized-energy mean-density-anomaly `G` modes with finite/zero/inactive endpoint limits, signed-unit normalization and signatures, exact constant-null handling, surface-referenced spectral `F` integration, aligned sampled transforms with direct `G` projection and diagnostic `F` synthesis, and rejection of numerically finite representations of infinite generalized eigenmodes

### Earlier V2 development

- added the parallel v2 `IMEigenvalueProblem` architecture with physical-coordinate EVP descriptors, structured operators, first-class boundary conditions, coordinate-aware spectral/WKB/density/finite-difference solvers, and `IMBasisSet`
- moved v2 modal normalization ownership to `IMBasisSet` so solved native modes can be reinterpreted with different normalizations without re-solving the EVP
- added v2 analytical solution families for constant and exponential stratification
- renamed v2 wave EVPs to `waveModesAtWavenumber` and `waveModesAtFrequency`, with `hydrostaticGModes` now using the fixed-frequency wave EVP at `omega=0`
- added native-coordinate v2 inner products for spectral, WKB spectral, and density spectral solvers, including endpoint surface and bottom weights outside the interior integral
- added separate `UnitTestsV2` regression coverage and `ExamplesV2` scripts for the new IM stack
- added index policies for active-boundary mode counts, including the manuscript positive-boundary and negative-boundary PE conventions
- added v2 regression tests for coordinate pullbacks, index policies, finite-difference grid ownership, basis-owned normalization, partial-depth Gram matrices, and unsupported analytical-basis operations
- expanded the partial-depth observational projection example with configurable sampling scenarios, Gaussian ensemble coefficient realizations, physical G-mode potential-energy spectra, observed-space completeness diagnostics, and spectral-window plots
- added a hydrostatic quadrature accuracy comparison example that separates mode-shape, root, solver, retained-band, and positive-weight effects in Parseval diagnostics
- added geostrophic normalization support and diagnostics for the analytical constant and exponential stratification solvers, with smoke-test coverage
- updated the physical sampling quadrature optimization example to use geostrophic normalization

## [1.4.0] - 2026-05-08
- added component-role-aware `InternalModesBasis`, `InternalModesTransform`, and `InternalModesProjection` classes for vertical transforms, fixed-grid projections, observational projections, vertical spectra, and annotated persistence
- added vertical-transform examples for wave-vortex-model integration sketches and component-role-aware modal spectra
- added a direct `ClassAnnotations ^1.2.1` dependency for NetCDF-backed annotated persistence through `CAAnnotatedClass`

## [1.3.0] - 2026-05-06
- ported the spectral speedup refactor into the `@InternalModesSpectral` class-folder layout
- added `modesAtQuadraturePoints` for computing resolved modes on mode-adapted quadrature points
- refreshed spectral internal state when `nEVP` changes, including WKB stretched-grid derivative state
- fixed `BSpline` initialization in `InternalModesBase`
- added speedup exploration helpers under `Extras/SpeedupExploration`
- expanded automated smoke coverage for quadrature points, `nEVP` refresh behavior, legacy quadrature wrappers, and spline initialization

## [1.2.0] - 2026-04-20
- modernized the non-spectral constructors to explicit name-value initialization while preserving the `InternalModesSpectral` and `InternalModesWKBSpectral` constructor contracts
- updated the `InternalModes` wrapper to forward constructor options through the new named-argument paths without breaking wrapper-level properties such as diagnostics and boundary-condition settings
- renamed the former script-only `UnitTests` collection to `Examples`, repaired the example scripts and local asset loading, and added a real `matlab.unittest` smoke suite under `UnitTests`
- refreshed direct constructor call sites and documentation snippets to match the new constructor forms and current spline dependency APIs
- modernized the core user-facing mode API to lowerCamel method names such as `modesAtFrequency`, `modesAtWavenumber`, `surfaceModesAtWavenumber`, and `showLowestModesAtFrequency`, while preserving the legacy UpperCamel names as hidden compatibility aliases and updating the docs, examples, and smoke tests accordingly
- modernized the direct eigenmatrix helper API to lowerCamel `eigenmatricesForFrequency` and `eigenmatricesForWavenumber`, while preserving the legacy UpperCamel names as hidden compatibility aliases and updating the developer-facing docs and smoke tests accordingly

## [1.1.0] - 2026-04-09
- updated spline-based interpolant construction to the `Distributions` 2.0 named-argument API
- raised the `SplineCore` dependency floor to `^2.0` and added a direct `Distributions ^2.0` dependency for the package's direct `NormalDistribution` usage

## [1.0.1] - 2025-12-09
- Initial CI release
